from VCFclass import VCF
import argparse

parser = argparse.ArgumentParser(description='Apply within-position filter to VCF file of called mutations.')

parser.add_argument('--VCF', metavar='VCF file', type=str,
                    help='VCF file of mutation calls')
parser.add_argument('--BAMfiles', metavar='BAMfile of each sample', nargs='+',
                    help='''BAM files to annotate from. One BAM file per sample must be specified. BAM files may be entered in one of two ways:
                            A single file (if one sample is contained in the VCF (e.g., --BAMfiles sample_1.bam)
                            A list of files in the same order as your VCF samples (e.g., --BAMfiles sample_1.bam sample_2.bam sample_3.bam)''')
parser.add_argument('--minimum_frequency', metavar='min freq', type=float, default=0.01,
                    help="minimum frequency of SNPs; SNPs below this frequency will be removed. Default is %(default)s (aka 1%%)")
parser.add_argument('--min_read_distance', metavar='Minimum alt location', type=float, default=0,
                    help='in-read position cutoff. If average alt location is further than X from the end of the read, variant is filtered out. X is a float value from 0 to 1. Default %(default)s (i.e., no read position filter)')
parser.add_argument('--min_base_quality', metavar='Minimum Base Quality', type=int, default=30,
                    help='Minimum base quality of reads to be considered for this filter. Default %(default)s.')
parser.add_argument('--min_mapping_quality', metavar='Minimum Mapping Quality', type=int, default=10,
                    help='Minimum mapping quality of reads to be considered for this filter. Default %(default)s.')
parser.add_argument('--logfile', metavar='Log file save location', type=str, default='snp_filter_log.log',
                    help='Location to save log file. Default %(default)s')
parser.add_argument('--stat_test',help='Which test to use to determine if alt bases are in a significantly different position than ref bases. Default %(default)s.', 
                    choices=['t-test', 'mann-whitney', 'bootstrap'], default='bootstrap')
parser.add_argument('--no_frequency_adjustment', help='Don\'t adjust frequencies of SNPs that are kept', action='store_false')
parser.add_argument('--keep_filtered_snps', help='Don\t remove SNPs that fail filter; instead, mark as a filter FAIL', action='store_true')

args = parser.parse_args()

print('Loading VCF file...')
vcf = VCF(args.VCF, bamfiles=args.BAMfiles)

print('Beginning filtering process...')
vcf.apply_position_filter(removeFails=args.keep_filtered_snps, freq_cutoff=args.minimum_frequency, 
                          in_read_cutoff=args.min_read_distance, min_base_quality=args.min_base_quality,
                          min_mapping_quality=args.min_mapping_quality, log_log=args.logfile,
                          adjust_freq=args.no_frequency_adjustment, stats_method=args.stat_test)

print(f'Saving VCF to {args.VCF[:-4]+"_filtered.vcf"}')
vcf.to_vcf(args.VCF[:-4]+'_filtered.vcf')
