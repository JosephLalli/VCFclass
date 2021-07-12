.. highlight:: python

=========================================
Typical Workflow
=========================================

Opening a file
==============

To begin, import the VCFclass module and open a
:class:`VCFclass.AlignmentFile`::

   from VCFclass import VCF
   vcf_file = VCF('myvcf.vcf')

The above command opens the file :file:`myvcf.vcf` and saves
it to a variable called vcf_file. VCFclass will open :file:`myvcf.vcf`,
detect the VCF's header, samples, and mutation calls, and load these
data into memory. This can take some time for large VCF files.

Depending on your downstream applications, it can be useful to include
a reference fasta, a gtf file location, or the location of the 
sam/bamfile used to create the vcf file. These can be optionally specified 
in the inital VCF import::

   myVCF = VCF('vcf_file.vcf', 
               refFile='my_ref.fasta',
               gtfFile='my_gtf.gtf',
               bamfiles='my_bamfile.bam')


Annotating a VCF file
=====================

To assign each mutation a gene and synonymous/non-synonymous status from a GTF annotation file and reference Fasta, simply run the following:
::

   myVCF = myVCF.annotate()
   
If you did not specify either a GTF file or a fasta file when you initially created the VCF object, you can pass these GTF annotation file as arguments at the time of annotation.
::

   myVCF = myVCF.annotate(gtf_file='my_gtf.gtf', ref_file='my_ref.fasta')

Assign BAM files to VCF file
===========================
Some downstream functions (such as in-read position filtering) require access to a sample's BAM file.
Bamfiles can be provided through :meth:`VCF.add_bamfile_locations` in one of three ways.

If only one sample is present in the VCF::

   myVCF.add_bamfile_locations('sample_bamfile.bam')

If multiple samples are present in the VCF, and you are sure of their
order in the VCF file, you can pass a list of BAM filenames. This list must be in the same order as the samples in your VCF file.
::

   myVCF.add_bamfile_locations(['sample1_bamfile.bam', 'sample2_bamfile.bam'])

If multiple samples are present in the VCF and you are *not* sure of their
order in the VCF file, you can pass a dictionary of the format {sample_name: bamfile}:
::

   myVCF.add_bamfile_locations({'sample1': 'sample1_bamfile.bam', 'sample2': 'sample2_bamfile.bam'})

Alternately, the string, list, and dictionary above can be passed to the VCF 

Applying the in-read position filter to a VCF file
==================================================

First, annotate the VCF and add to the VCF. Then run :meth:`VCFclass.VCF.apply_position_filter`::

   myVCF.apply_position_filter()

The default options should apply most of the time. For more information see :meth:`VCFclass.VCF.apply_position_filter`.


Combining multiple VCF files into one multi-sample VCF
======================================================

::

   myVCF.merge(my_other_VCF)

Averaging VCF files from technical replicates
=============================================

::

   myVCF.average(my_other_VCF)

Exporting VCF data as Pandas dataframe
======================================

::

   myVCF_DF = myVCF.to_dataframe()

Export per-nucleotide read counts to numpy array
================================================

::

   myVCF_array = myVCF.to_numpy()
