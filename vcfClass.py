import pandas as pd
import os,sys
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import pysam
from collections import Counter, defaultdict
from operator import itemgetter
from scipy import stats
import numpy as np


def importVCF(location):
	return VCF(location)

class VCF:
	def __init__ (self, location, refFile=None, gtfFile=None, bamfiles={}):
		self.vcfFileName = location

		#handle gz files:
		try:
			with open(location, 'r') as vcffile:#imgigi
				self._vcflines = vcffile.readlines()
		except:
			try:
				with gzip.open(location, 'r') as vcffile:#imgigi
					self._vcflines = vcffile.readlines()
			except:
				raise

		self._rowstoskip = self._getVCFStart()

		self.header = Header(self._vcflines[0:self._rowstoskip-1])
		self.reference = {}
		if refFile:
			self.addReferenceFile(refFile)
		self.gtffile = gtfFile
		self.refFile = refFile
		self.bamfiles={}

		self.samples = []
		self.samples = self._vcflines[self._rowstoskip-1].strip().split('\t')[9:]

		self.mutations = [MutCall(row, [sample for sample in self.samples]) for row in self._vcflines[self._rowstoskip:]]
		self.SNPs = [mut for mut in self.mutations if mut.type == 'SNP']
		self.indels = [mut for mut in self.mutations if mut.type == 'insertion' or mut.type == 'deletion']
		self._hashedmuts = {mut.chrom:{int(mut1.pos):mut1 for mut1 in self.mutations if mut1.chrom == mut.chrom} for mut in self.mutations}
		self.read_tsv_args = {'sep':'\t', 'keep_default_na':False, 'na_values':['-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN', '#N/A','N/A', '#NA', 'NULL', 'NaN', '-NaN', 'nan', '-nan']}

		self.add_bamfile_locations(bamfiles)

	def _getVCFStart(self):
		rowstoskip = 0
		#Opens stresentative as a text file and identifies the row in which our data begins
		for index, line in enumerate(self._vcflines):
			if "#CHROM" in line:
				rowstoskip = index+1
				break
		return rowstoskip

	def addReferenceFile(self, ref):
		self.refFile = ref
		with open(ref, 'r') as r:
			ref = SeqIO.parse(r, 'fasta')
			self.reference = {seq.id:str(seq.seq) for seq in ref}

	def averageWithVCF(self,otherVCF,newSampleName=None):
		'''
		*assumes one sample per vcf*
		take intersection of mutations and average their coverages/frequencies etc.
		'''
		print (f'Averaging {self} with {otherVCF}')
		
		mutsnotinotherVCF = {(mut.chrom,mut.pos) for mut in self.mutations} #set
		
		for mut in otherVCF:
			try:
				self._hashedmuts[mut.chrom][mut.pos] = self._hashedmuts[mut.chrom][mut.pos].average(mut)
				mutsnotinotherVCF.discard((mut.chrom,mut.pos))
				if self._hashedmuts[mut.chrom][mut.pos] == None:
					mutsnotinotherVCF.add((mut.chrom,mut.pos))
			except KeyError:
				pass
			
		for chrom, pos in mutsnotinotherVCF:
			self.removemut(chrom, pos)
		
		self.renameSample(self.samples[0],newSampleName)
		#self.header = self.header.combineHeaders(otherVCF.header)
		return self

	def averageSamples(self,sampleA,sampleB,newSampleName=None):
		'''take intersection of mutations and average their coverages/frequencies etc.'''

		intersection = [mut for mut in self.mutations if mut.hasSample(sampleA) and mut.hasSample(sampleB)]
		#iterate through all mutations
		for mut in intersection:
			mut = mut.averageSamples(sampleA,sampleB)
		
		self.renameSample(sampleA, newSampleName)
		self.deleteSample(sampleB)
		return self

	def renameSample(self, origname, newName):
		'''update name of sample origname in both list of samples (self.samples)
		and in all mutations'''
		if newName == None:
			return 1
		for mut in self.mutations:
			mut.renameSample(origname, newName)
		self.samples = [sample if sample != origname else newName for sample in self.samples]

	def deleteSample(self, samplename):
		for mut in self.mutations:
			mut.deleteSample(samplename)
		self.samples.remove(samplename)
	
	def mergeVCFs(self,otherVCF):
		for mut in otherVCF:
			try:
				self._hashedmuts[mut.chrom][mut.pos].addSamples(mut.samples, mut)
				self.samples.extend(mut.samples)
			except KeyError:
				self.addMut(mut)
		return self

	def addMut(self,newMut):
		raise Exception('addMut not implemented yet')

	def removemut(self, chrom, pos):
		'''muts are stored in self.mutations and self._hashedmuts. 
		this removes all muts w/ chrom and pos from both lists.'''
		try:
			print('deleting ' + chrom + ' '+str(pos))
			self.mutations = [mut for mut in self.mutations if ((mut.chrom != chrom) or (mut.pos != pos))]
			del self._hashedmuts[chrom][pos]
		except KeyError:
			pass

	def annotate(self, gtffile=None, reference=None):
		if not self.gtffile:
			self.gtffile = gtffile
		if not self.refFile:
			self.addReferenceFile(reference)

		coding_regions = extractCodingRegions(gtffile)

		# use gene coordinates to create coding sequences from reference sequences
		transcripts = createTranscripts(coding_regions, self.reference)
		AAtranscripts = {gene:str(Seq(transcript).translate()) for gene, transcript in transcripts.items()}

		for segment in coding_regions.keys():
			for gene in coding_regions[segment].keys():
				priorExonLength = 0
				for start, stop in coding_regions[segment][gene]:
					for mut in self.fetchSNPs(segment, start, stop+1):
						offset = start - priorExonLength
						mut.annotate(gene, offset, transcripts[gene], AAtranscripts[gene])
					priorExonLength += (stop+1-start)

	def add_bamfile_locations(self, bamfiles):
		if len(bamfiles) == 0:
			self.bamfiles = {}
		elif type(bamfiles) == str and len(self.samples) == 1:
			bamfile = os.path.abspath(bamfiles)
			self.bamfiles = {self.samples[0]:bamfile}
		elif type(bamfiles) == list and len(bamfiles) == len(self.samples):
			self.bamfiles = {sample:os.path.abspath(bamfile) for sample, bamfile in zip(self.samples, bamfiles)}
		elif type(bamfiles == dict):
			for sample, bamfile in bamfiles.items():
				if sample not in self.samples:
					raise Exception(f'Sample {sample} not in VCF.')
			self.bamfiles.update(bamfiles)
		else:
			raise Exception ('Unable to associate bamfiles with samples.')

	def fetchSNPs(self, chrom, start=0, end=float('inf')):
		return [mutation for mutation in self.mutations \
			if (mutation.chrom == chrom and mutation.pos >= start and mutation.pos < end)]

	# def benjamini_hochberg_pfilter(self, pvalues, alpha=0.05):
		
	# def apply_position_filter_to_sample(self, sample)
	def apply_position_filter(self, removeFails = True, signifigance_at = 0.05, removeFailsMethod = 'Bonferroni', in_read_cutoff=0.1, freq_cutoff=0.01):
		#count mutations deleted for fun
		if removeFailsMethod == 'Bonferroni':
				#Number of comparisons that will be made is number of mutations * number of samples
				# cutoff = signifigance_at/(len(self.mutations)*len(self.samples))
				df = self.to_dataframe()
				self.pval_cutoff = signifigance_at/len(df.loc[df.FREQ > 0.01])
		else:
			self.pval_cutoff = signifigance_at

		failcount = 0
		# filter_sample = functools.partial(apply_position_filter_to_sample, removeFails = removeFails, cutoff=cutoff, in_read_cutoff=0.1)
		# with mp.Pool(8) as pool:
		# 	pool.map(apply_position_filter_to_sample, )
		for sample in self.samples:
			try:
				bamfile = self.bamfiles[sample]
			except KeyError:
				print (f'''Position filter requires the original bamfile from which SNPs were called.\n
						   Bamfile for sample {sample} is missing.)''')
				print ("You can add the bamfile by using vcf.add_bamfile_locations(format: {sample:bamfileLocation})")
			bam = pysam.AlignmentFile(bamfile, 'rb')
			print (f'Processing {sample}')
			origMuts = len(self.mutations)
			for mut in self.mutations:
				if mut.hasSample(sample):
					samp = mut.get(sample)
					RD = samp.RD
					AD = samp.AD
					if .5-abs(.5-samp.freq) > freq_cutoff:
						pos_filter_result = mut.apply_pos_filter(bam, sample, self.pval_cutoff, removeFails, in_read_cutoff)
					else:
						continue
					if pos_filter_result != 'PASS':
						print(pos_filter_result)
					if pos_filter_result == 'FAIL':
						failcount += 1
						print (failcount)
						if not mut.stillHasSNPs():
							self.removemut(mut.chrom, mut.pos)
		print (f'Finished.\n{failcount} SNVs failed position filter.')
		return self

	def get_synonymous(self):
		return [mutation for mutation in self.mutations \
			if mutation.AAtype == 'Synonymous']

	def get_nonsynonymous(self):
		return [mutation for mutation in self.mutations \
			if mutation.AAtype == 'Nonsynonymous']

	def to_dict(self):
		'''returns list of dictionaries. Each dictionary is 'property':value for each snp in each sample.'''
		return [entry for mut in self.mutations for entry in mut.to_dict()]

	def to_numpy(self, referenceFile=None, cutoff_freq=0):
		# def generate_ref_numpy(self, referenceFile):
		if len(self.reference) < 1:
			if referenceFile:
				self.addReferenceFile(referenceFile)
			else:
				print('This function requires VCF to have a reference sequence.')
				return None

		# self.refFile = ref
		# with open(ref, 'r') as r:
		# 	refseq = SeqIO.parse(r, 'fasta')
		# 	self.reference = {seq.id:str(seq.seq) for seq in ref}

		# refseq = list(SeqIO.parse(refseq, 'fasta'))

		concatrefseq = ""
		segStarts = dict()
		segCoords = dict()
		runningtally = 0
		for chrom, seq in self.reference.items():
			segStarts[chrom.split('_')[-1]] = runningtally
			segCoords[chrom.split('_')[-1]] = (runningtally, runningtally+len(seq))
			runningtally += len(seq)
			concatrefseq += seq

		totalConcatLength = len(concatrefseq)

		df = self.to_dataframe()

		#only report muts that are within cutoff_freq; adjust all other read counts to 0/1
		#first adjust SNPs under cutoff_freq (all ref)
		df.loc[(df.FREQ > 0) & (df.FREQ < cutoff_freq), 'AD'] = 0
		df.loc[(df.FREQ > 0) & (df.FREQ < cutoff_freq), 'DP'] = df.loc[(df.FREQ > 0) & (df.FREQ < cutoff_freq), 'RD']
		df.loc[(df.FREQ > 0) & (df.FREQ < cutoff_freq), 'FREQ'] = 0
		# then adjust SNPs above 1-cutoff_freq (all alt)
		df.loc[(df.FREQ < 1) & (df.FREQ > 1-cutoff_freq), 'RD'] = 0
		df.loc[(df.FREQ < 1) & (df.FREQ > 1-cutoff_freq), 'DP'] = df.loc[(df.FREQ < 1) & (df.FREQ > 1-cutoff_freq), 'AD']
		df.loc[(df.FREQ < 1) & (df.FREQ > 1-cutoff_freq), 'FREQ'] = 1

		df.chrom = df.chrom.str.split('_').str[-1]

		df['inConcatPos'] = df.pos
		for seg, offset in segStarts.items():
			df.loc[df.chrom == seg,'inConcatPos'] += offset

		nucDict={'A':0,'C':1,'G':2,'T':3}
		df['ref_nuc'] = df.ref.map(nucDict)
		df['alt_nuc'] = df.alt.map(nucDict)

		RDdf = df[['sampleID','inConcatPos','ref_nuc','RD', 'chrom']]
		RDdf = RDdf.sort_values('RD',ascending=False).reset_index(drop=True)
		RDdf = RDdf.loc[~RDdf[['sampleID','inConcatPos','chrom','ref_nuc']].duplicated()]
		readCtsDF = pd.pivot_table(RDdf[['sampleID','inConcatPos','ref_nuc','RD']], columns='ref_nuc', values = 'RD', index=('sampleID', 'inConcatPos'))

		ADdf = df[['sampleID','inConcatPos','alt_nuc','AD', 'chrom']]
		ADdf = ADdf.sort_values('AD',ascending=False).reset_index(drop=True)
		ADdf = ADdf.loc[~ADdf[['sampleID','inConcatPos','chrom','alt_nuc']].duplicated()]
		readCtsDF.update(pd.pivot_table(df[['sampleID','inConcatPos','alt_nuc','AD']], columns='alt_nuc', values = 'AD', index=('sampleID', 'inConcatPos')))

		readCtsDF = readCtsDF.unstack().fillna(0)
		positions = readCtsDF[0].columns
		readCts = np.zeros((readCtsDF.shape[0],totalConcatLength,4))
		readCts[:,positions,0] = readCtsDF[0].to_numpy()
		readCts[:,positions,1] = readCtsDF[1].to_numpy()
		readCts[:,positions,2] = readCtsDF[2].to_numpy()
		readCts[:,positions,3] = readCtsDF[3].to_numpy()

		samplelist = list(readCtsDF.index.get_level_values(0))

		return readCts, samplelist, list(positions)

	def to_dataframe(self):
		'''export VCF as tidy dataframe'''
		return pd.DataFrame(self.to_dict()).rename(columns={'sample':'sampleID'}) #sample is a reserved term in Pandas

	def to_vcf(self, location):
		with open (location, 'w') as outfile:
			outfile.write(str(self.header))
			columnheader = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+'\t'.join(self.samples)+'\n'
			outfile.write(columnheader)
			outfile.writelines([str(mutation)+'\n' for mutation in self.mutations])

	def __len__(self):
		return len(self.mutations)

	def __iter__(self):
		for mutation in self.mutations:
			yield mutation

	def __str__(self):
		return f"VCF containing {len (self.samples)} samples and {len(self.mutations)} mutation calls"


class MutCall:
	def __init__(self, row, samples):
		self._rawlist = row.split('\t')
		self.samples = samples
		self.chrom = self._rawlist[0]
		#internal representation will be 0-indexed. Will export VCFs as 1-indexed, all others exported as 0-indexed.
		self.pos = int(self._rawlist[1])-1
		self.id = self._rawlist[2]
		self.ref = self._rawlist[3]
		self.alt = self._rawlist[4]
		self.qual = self._rawlist[5]
		self.filter = self._rawlist[6]
		self.info  = self._rawlist[7]
		self.format = self._rawlist[8]
		self.type = self._determinetype()
		self._sampledata = {sample:SampleMut(sample, data, self.format) for data, sample in zip(self._rawlist[9:], self.samples)}

		del self._rawlist

		self.genes = []
		self.genePOS = []
		self.AArefs = []
		self.AAalts = []
		self.AApos = []
		self.AAstrs = []
		self.AAtypes = []

	def _determinetype(self):
		if len(self.ref) > 1:
			return "deletion"
		elif len(self.alt) > 1:
			return "insertion"
		elif len(self.alt) == 1:
			return "SNP"
		else:
			return "translocation"

	def get(self, sample):
		return self._sampledata[sample]

	def stillHasSNPs(self):
		return np.array([self.hasSample(sample) for sample in self.samples]).any()

	def bootstrap(self, array, bootsize):
		return np.random.choice(array, len(array)*bootsize, replace=True).reshape(-1, bootsize).mean(axis=0)

	def apply_pos_filter(self, bam, sample, cutoff=0.05, removeFails=True, in_read_cutoff=0.1, min_base_quality=30, min_mapping_quality=10, log_loc='snp_filter_log.log'):
		'''given pysam alignment object and relevent sample, determine whether average
		in-read positions for ref and alt alleles come from different distributions.
		If so, it sets that sample's minority allele read ct to 0 (since these are presumed to be alignment errors).
		If that means the ref allele is 100%, the sample is removed from the mutation.
		If that means the alt allele is 100%, the sample is retained with updated statistics.'''

		chrom = self.chrom
		pos = self.pos
		pileup = bam.pileup(contig=chrom, start=pos, end=pos+1, truncate=True, stepper="nofilter")
		column = next(pileup)
		sampmut = self.get(sample)

		column.set_min_base_quality(min_base_quality)
		bases = [base.upper() for base in column.get_query_sequences()]
		map_quals = column.get_mapping_qualities()

		soft_clipped_positions = [read.alignment.get_reference_positions().index(pos) if pos in read.alignment.get_reference_positions() else 0 for read in column.pileups]
		read_lengths = [read.alignment.query_alignment_length for read in column.pileups]

		positions = defaultdict(list)
		for position, read_length, base, mQ in zip(soft_clipped_positions, read_lengths, bases, map_quals):
			if mQ < min_mapping_quality:
				continue
			try:
				inReadPos = position/read_length
				positions[base].append(.50-abs(inReadPos-.50)) #distance of SNP from end of read
			except:
				print (base)
				print (bases)

		avg_ref_pos = np.mean(positions[self.ref.upper()])
		avg_alt_pos = np.mean(positions[self.alt.upper()])

		if (self.ref.upper() not in positions.keys()) or (self.alt.upper() not in positions.keys()):
			p_value = 1 #Definitely not an error. p-value = 1.
		elif (len(positions[self.ref.upper()]) <= 1) or (len(positions[self.alt.upper()]) <= 1):
			p_value = 1 #Can't do pvalue calc on one position
		else:
			# ref_positions = np.array(positions[self.ref.upper()])
			# alt_positions = np.array(positions[self.alt.upper()])
			# ref_positions = np.unique(np.array(positions[self.ref.upper()]))
			# alt_positions = np.unique(np.array(positions[self.alt.upper()]))
			ref_positions = self.bootstrap(np.array(positions[self.ref.upper()]), int(1/cutoff)+1)#np.array([np.mean(np.random.choice(np.array(positions[self.ref.upper()]), len(positions[self.ref.upper()]), replace=True)) for _ in range(int(1/cutoff)+1)])
			alt_positions = self.bootstrap(np.array(positions[self.alt.upper()]), int(1/cutoff)+1)#np.array([np.mean(np.random.choice(np.array(positions[self.alt.upper()]), len(positions[self.alt.upper()]), replace=True)) for _ in range(int(1/cutoff)+1)])
			# ref_positions = np.random.choice(np.array(positions[self.ref.upper()]), 5000, replace=True)
			# alt_positions = np.random.choice(np.array(positions[self.alt.upper()]), 5000, replace=True)
			bigger_than_pvalue = np.mean(ref_positions>alt_positions)
			less_than_pvalue = np.mean(ref_positions<alt_positions)
			p_value = min(bigger_than_pvalue, less_than_pvalue)
			# p_value = np.mean(np.abs(ref_positions-alt_positions) < 0.05)
			# print(p_value)
			# tstat, p_value = stats.mannwhitneyu(ref_positions, alt_positions)

			# tstat, p_value = stats.mannwhitneyu(upsampled_ref, upsampled_alt)
			# tstat, p_value = stats.ttest_ind(positions[self.ref.upper()], positions[self.alt.upper()],equal_var=False)

		sampmut.avg_ref_pos = avg_ref_pos
		sampmut.avg_alt_pos = avg_alt_pos
		sampmut.read_pos_p_value = p_value

		# regardless of "ref" and "alt", we want to examine the in-read position of the allele w/
		# fewer mapped entries.
		if len(positions[self.ref.upper()]) < len(positions[self.alt.upper()]):
			dominant_allele_positions = positions[self.alt.upper()]
			minor_allele_positions = positions[self.ref.upper()]
		else:
			minor_allele_positions = positions[self.alt.upper()]
			dominant_allele_positions = positions[self.ref.upper()]

		location_cutoff = in_read_cutoff

		# Read position logic:
		# If the positions are significantly different and are separated by a distinct amount, the mutation is not valid.
		# If the positions are not significantly different but they are on average
		# within last 10% of end of read, the mutation is not valid.
		close_to_end_of_read = np.mean(minor_allele_positions) < location_cutoff
		sig_diff_pos_from_major_allele = (p_value < cutoff)
		# print()
		actually_different_positions = True #np.abs(np.mean(minor_allele_positions)-np.mean(dominant_allele_positions)) > 0.05

		if sig_diff_pos_from_major_allele and actually_different_positions:
			sampmut.position_filter = 'FAIL'
			with open(log_loc, 'a') as log:
				print (f'For sample {sample}, mutation {self.chrom.split("_")[-1]} {self.pos} failed read position filter.',file=log)
				print (f'Avg ref position: {np.round(avg_ref_pos,3)}(n={len(positions[self.ref.upper()])}). Avg alt position: {np.round(avg_alt_pos,3)}(n={len(positions[self.alt.upper()])}).',file=log)
				print (f'Read location cutoff was {np.round(location_cutoff, 4)}. p-value {p_value} is less than cutoff {cutoff}.', file = log)
			#Write to screen
			freq= len(positions[self.ref.upper()])/(len(positions[self.ref.upper()])+len(positions[self.alt.upper()]))
			if (freq > 0.01) and (freq < 0.99):
				print (f'For sample {sample}, mutation {self.chrom.split("_")[-1]} {self.pos} read position filter: {sampmut.position_filter}.')
				print (f'Avg ref position: {np.round(avg_ref_pos,3)}(n={len(positions[self.ref.upper()])}). Avg alt position: {np.round(avg_alt_pos,3)}(n={len(positions[self.alt.upper()])}).')
				print (f'Read location cutoff was {np.round(location_cutoff, 4)}. p-value {p_value} is less than cutoff {cutoff}.')

			if removeFails:
				if sampmut.RD>sampmut.AD:
					sampmut.zeroOut()
				elif sampmut.RD<=sampmut.AD:
					sampmut.position_filter = 'PASS'
					sampmut.DP -= sampmut.RD
					sampmut.RD = 0
					sampmut.freq = 1.0
					sampmut.freqstring='100%'
					sampmut.PVAL = np.round(p_value,4)
					sampmut.RBQ = 0
					sampmut.RDF = 0
					sampmut.RDR = 0
					sampmut.update()
			print (sampmut.position_filter, sampmut.freqstring)
		else:
			sampmut.position_filter = 'PASS'

		return sampmut.position_filter

	def averageSamples(self, sampleA, sampleB):
		self._sampledata[sampleA] = self._sampledata[sampleA].average(self._sampledata[sampleB])
		return self
	
	def renameSample(self,oldsamplename,newName):
		self._sampledata[newName] = self._sampledata.pop(oldsamplename)
		self.samples = [sample if sample != oldsamplename else newName for sample in self.samples]

	def deleteSample(self, sampletoDelete):
		self._sampledata.pop(sampletoDelete)
		self.samples.remove(sampletoDelete)

	def average(self, otherMut):
		sample = self.get(self.samples[0])
		sample2 = otherMut.get(otherMut.samples[0])
		if sample2.SDP == 0 or sample.SDP==0:
			return None
		else:
			sample = sample.average(sample2)

		tempinfo = self.info.split(';')
		for ADPloc, item in enumerate(tempinfo):
			if "ADP" in item:
				ADP = int(item.split('=')[-1])
				for otheritem in otherMut.info.split(';'):
					if "ADP" in otheritem:
						otherADP = int(otheritem.split('=')[-1])
						break
				break
		
		tempinfo[ADPloc] = str(int((ADP+otherADP)/2))
		self.info = "ADP="+";".join(tempinfo)
		return self
	
	def hasSample(self, samplename):
		return self._sampledata[samplename].exists()
	
	def addSamples(self,samplelist,newMut):
		pass

	def annotate(self, gene, adjust, dnaseq, AAseq):
		inGenePos = self.pos-adjust #pos is 0-indexed, so this is 0-indexed.
		self.genes.append(gene)
		self.genePOS.append(inGenePos)
		if self.ref != dnaseq[inGenePos]:
			print (self.pos)
			print (self.ref)
			print (inGenePos)
			print (adjust)
			print (gene)
			print (dnaseq)
			print (AAseq)
			print (dnaseq[inGenePos])
			raise Exception
		
		#Calc assuming type is SNP:
		AApos = int((inGenePos)/3)
		refAA = AAseq[AApos]
		inAApos = (inGenePos)%3
		refcodon = dnaseq[inGenePos-inAApos:inGenePos-inAApos+3]
		altcodon = refcodon[:inAApos]+self.alt+refcodon[inAApos+1:]
		altAA = str(Seq(altcodon).translate())

		#currently not calling indels, this is incomplete (cannot handle partial indels, multiple )
		if self.type == 'deletion':
			refAA = '-'*int(len(self.ref)/3)
			if len(refAA) == 0:
				refAA = 'd'
		if self.type == 'insertion':
			altAA = '-'*int(len(self.ref)/3)
			if len(refAA) == 0:
				altAA = 'd'

		self.AApos.append(AApos)
		self.AArefs.append(refAA)
		self.AAalts.append(altAA)
		self.AAstrs.append(refAA+str(AApos)+altAA)
		if refAA == altAA:
			self.AAtypes.append('Synonymous')
		else:
			self.AAtypes.append('Nonsynonymous')

	def to_dict(self):
		selfDict = []
		for samplename in self.samples:
			sample = self.get(samplename)
			result = {'sample':samplename,'chrom':self.chrom, 'pos':self.pos, 'id':self.id, 'ref':self.ref,'alt':self.alt,'qual':self.qual}
			result.update(sample.to_dict())
			if len(self.AApos) > 0: #if annotated
				for gene, inGenePos, AAref, AApos, AAalt, AAstr, AAtype in zip(self.genes, self.genePOS, self.AArefs, self.AApos,self.AAalts,self.AAstrs, self.AAtypes):
					geneSpecificResult = result.copy()
					geneSpecificResult.update({'gene':gene,'inGenePos':inGenePos,'refAA':AAref, 'AApos':AApos, 'altAA':AAalt, 'AAstr':AAstr, 'AAtype':AAtype})
					selfDict.append(geneSpecificResult)
			else:
				selfDict.append(result)
		return selfDict

	def __str__(self):
		#Add 1 to pos to account for 1-indexing of VCF files
		return '\t'.join([self.chrom, str(self.pos+1), self.id, self.ref,self.alt,self.qual,self.filter,self.info,self.format])+'\t'+'\t'.join(str(self.get(sample)) for sample in self.samples)
		
	def __iter__(self):
		for sample in self.samples:
			yield self._sampledata[sample]

class SampleMut:
	def __init__(self, samplename, data, formatinfo):
		self.name = samplename
		self.format = formatinfo
		self.other=[]
		
		for label, item in zip(self.format.split(":"), data.split(":")):
			item = item.strip('\n').strip('\"')
			
			if item == '.':
				item = 0
			elif item.isnumeric():
				item = int(item)
				
			if label == "GT":
				self.GT = item
			elif label == "GQ":
				self.GQ = item
			elif label == "SDP":
				self.SDP = item
			elif label == "DP":
				self.DP = item
			elif label == "RD":
				self.RD = item
			elif label == "AD":
				self.AD = item
			elif label == "FREQ":
				self.freqstring = item
				try:
					self.freq = round(float(item.rstrip('%'))/100,4)
				except:
					self.freq = round(float(self.freqstring)/100,4)
			elif label == "PVAL":
				self.PVAL = item
			elif label == "RBQ":
				self.RBQ = item
			elif label == "ABQ":
				self.ABQ = item
			elif label == "RDF":
				self.RDF = item
			elif label == "RDR":
				self.RDR = item
			elif label == "ADF":
				self.ADF = item
			elif label == "ADR":
				self.ADR = item
			else:
				self.other.append((label,item))

		
		# #When calculating freq and depth, Varscan removes AD reads due to quality filter, but not RD reads.
		# #But it reports only quality RD reads.
		# #That's a problem when you've got 7 crummy RD reads and a 100 good AD reads!
		# #I'll recalc depth and frequency here.
		if self.SDP != 0:
			self.DP = self.RD + self.AD
			self.freq = round(self.AD/self.DP, 4)
			self.freqstring = str(round((self.freq*100),4))+'%'

		self._properties = [self.GT,self.GQ,self.SDP,self.DP,self.RD,self.AD,self.freqstring,self.PVAL,self.RBQ,self.ABQ,self.RDF,self.RDR,self.ADF,self.ADR]

	def update(self):
		self._properties = [self.GT,self.GQ,self.SDP,self.DP,self.RD,self.AD,self.freqstring,self.PVAL,self.RBQ,self.ABQ,self.RDF,self.RDR,self.ADF,self.ADR]

	def average(self, otherSample):
		self.freq = round((self.freq+otherSample.freq)/2, 4)
		self.freqstring = str(self.freq*100)+"%"
		self.GQ = int(round((self.GQ+otherSample.GQ)/2, 0))
		self.SDP += otherSample.SDP
		self.DP += otherSample.DP
		oldAD = self.AD
		self.AD = int(round(self.DP*self.freq,0))
		oldRD = self.RD
		self.RD = self.DP-self.AD
		self.RBQ = int(round((self.RBQ+otherSample.RBQ)/2, 0))
		self.ABQ = int(round((self.ABQ+otherSample.ABQ)/2, 0))

		try:
			RDFtoRD = ((self.RDF/oldRD)+(otherSample.RDF/otherSample.RD))/2
			self.RDF = int(round(self.RD*RDFtoRD,0))
			self.RDR = self.RD-self.RDF
		except ZeroDivisionError:
			self.RDF=0
			self.RDR=0

		try:
			ADFtoRD = ((self.ADF/oldAD)+(otherSample.ADF/otherSample.AD))/2
			self.ADF = int(round(self.AD*ADFtoRD,0))
			self.ADR = self.AD-self.ADF
		except ZeroDivisionError:
			self.RDF=0
			self.RDR=0
		
		self.update()
		
		return self

	def exists(self):
		if self.SDP == 0:
			return False
		else:
			return True

	def zeroOut(self):
		'''If a sample fails a filter or otherwise needs to be deleted in only one sample,
		   it's zero-ed out so that it's blank. Importantly, deleting a sample from one mutation
		   is different than deleting a sample from the whole vcf.'''
		self.GT='./.'
		self.GQ=0
		self.SDP=0
		self.DP=0
		self.RD=0
		self.AD=0
		self.freqstring=0
		self.PVAL=0
		self.RBQ=0
		self.ABQ=0
		self.RDF=0
		self.RDR=0
		self.ADF=0
		self.ADR=0
		self.update()

	def to_dict(self):
		self.update()
		return {"GT":self.GT,"GQ":self.GQ,"SDP":self.SDP,"DP":self.DP,"RD":self.RD,"AD":self.AD,"FREQ":self.freq,"PVAL":self.PVAL,"RBQ":self.RBQ,"ABQ":self.ABQ,"RDF":self.RDF,"RDR":self.RDR,"ADF":self.ADF,"ADR":self.ADR}
	
	def __str__(self):
		self.update()
		return ":".join([str(item) for item in self._properties])

class Header:
	def __init__(self, headertext):
		self.text = headertext

	def combineHeaders(self,otherheader):
		for line in otherheader.text:
			if line not in self.text:
				self.text.append(line)
	
	def __str__(self):
		return "".join(self.text)


def extractCodingRegions(gtffile):
	with open(gtffile, 'r') as g:
		gtf = g.readlines()

	coding_regions = {} 
	for line in gtf:
		line = line.replace("/", "_")
		lineitems = line.split("\t")
		segment_name = lineitems[0]
		annotation_type = lineitems[2]
		start = int(lineitems[3]) - 1  # adding the -1 here for 0 indexing
		stop = int(lineitems[4]) - 1	# adding the -1 here for 0 indexing
		gene_name = lineitems[8]
		gene_name = gene_name.split(";")[0]
		gene_name = gene_name.replace("gene_id ","")
		gene_name = gene_name.replace("\"","")

		if annotation_type.lower() == "cds":
			if segment_name not in coding_regions:
				coding_regions[segment_name] = {}
				coding_regions[segment_name][gene_name] = [[start, stop]]
			elif segment_name in coding_regions and gene_name not in coding_regions[segment_name]:
				coding_regions[segment_name][gene_name] = [[start, stop]]
			elif gene_name in coding_regions[segment_name]:
				coding_regions[segment_name][gene_name].append([start, stop])
	
	return coding_regions

def createTranscripts(coding_regions, ref_segments):
	transcripts = {}
	for segment in coding_regions:
		for gene in coding_regions[segment]:
			transcripts[gene] = ""
			coordinates = coding_regions[segment][gene]  # define the coding regions for each gene
			for start, stop in coordinates:   # loop through start/stop sites in coding regions
				sequence_chunk = ref_segments[segment][start:stop+1]
				transcripts[gene] = transcripts[gene] + sequence_chunk	 # append each piece of the transcript together 
	return transcripts