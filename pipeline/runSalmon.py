#!/usr/bin/env python3

import argparse
import shutil
import subprocess
import sys

def salmon_in_path():
    return shutil.which('salmon') is not None

def salmon_version():
    p = subprocess.Popen(['salmon', '--version'], shell=False,
                         stderr=subprocess.PIPE, encoding='utf-8')
    return 'salmon ' + p.communicate()[1].strip()

def print_cmd(args):
    print(salmon_version())
    print('Command line used:', file=sys.stderr)
    print('\t' + ' '.join(args), file=sys.stderr)

def salmon_index(transcripts, index, kmer_length=31, perfect_hash=False,
                 sample_interval=1, threads=1):
    args = [
        'salmon', 'index',
        '-t', transcripts,
        '-i', index,
        '-p', str(threads),
        '--type', 'quasi',
        '-s', str(sample_interval)
    ]
    print_cmd(args)
    p = subprocess.Popen(args, shell=False)
    p.wait()

def salmon_quant(index, read1, read2=None, output=None, seq_bias=False, gc_bias=False,
                 threads=1, libtype=True):
    args = ['salmon', 'quant',
            '-l', libtype,
            '-i', index]
    if read2 is not None:
        args += ['-1'] + read1 + ['-2'] + read2
    else:
        args += ['-r'] + read1

    args += ['-o', output,
             '-p', str(threads)]

    print_cmd(args)
    p = subprocess.Popen(args, shell=False)
    p.wait()

def parse_args():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(dest='subcommand')
    subparsers.required = True

    # index
    index_parser = subparsers.add_parser('index',
                                         help='create a '
                                         'salmon index (only '
                                         'supports salmon\'s '
                                         'quasi-alignment)')
    index_parser.add_argument('-t', '--transcripts',
                              help='transcript fasta file',
                              metavar='FILE',
                              required=True)
    index_parser.add_argument('-k', '--kmer-length',
                              help='k-mer length that should be used for the '
                              'quasi index (default: 31)',
                              type=int,
                              metavar='k',
                              default=31)
    index_parser.add_argument('-i', '--index',
                              help='salmon index',
                              metavar='DIR',
                              required=True)
    index_parser.add_argument('--perfect-hash',
                              help='build the index using a perfect hash '
                              'rather than a dense hash',
                              action='store_true')
    index_parser.add_argument('-s', '--sample-interval',
                              help='the interval at which the suffix array '
                              'should be sampled (default: 1)',
                              type=int,
                              default=1)
    index_parser.add_argument('-p', '--threads',
                              type=int, default=1,
                              metavar='N',
                              help='number of threads to use (default: 1)')


    # quant
    quant_parser = subparsers.add_parser('quant',
                                         help='quantify transcription '
                                         '(only supports salmon\'s '
                                         'quasi-alignment)')
    quant_parser.add_argument('-i', '--index',
                              help='salmon index',
                              metavar='DIR',
                              required=True)
    quant_parser.add_argument('-1', '--mates1',
                              help='file(s) containing the #1 mates',
                              metavar='FILE', nargs='+')
    quant_parser.add_argument('-2', '--mates2',
                              help='file(s) containing the #2 mates',
                              metavar='FILE', nargs='+')
    quant_parser.add_argument('-r', '--unmated-reads',
                              help='list of files containing single end reads',
                              metavar='FILE', nargs='+')
    quant_parser.add_argument('-o', '--output',
                              help='output directory',
                              metavar='DIR',
                              required=True)
    quant_parser.add_argument('--seq-bias',
                              help='perform sequence-specific bias correction',
                              action='store_true')
    quant_parser.add_argument('--gc-bias',
                              help='perform fragment GC bias correction',
                              action='store_true')
    quant_parser.add_argument('-l', '--libtype',
                              help='library type; see salmon documentation '
                              '(default: A)',
                              default='A')
    quant_parser.add_argument('-p', '--threads',
                              type=int, default=1,
                              metavar='N',
                              help='number of threads to use (default: 1)')

    args = parser.parse_args()

    if args.unmated_reads is not None and \
            ( args.mates1 is not None or args.mates2 is not None ):
        parser.error('either use pair end reads or single end reads, not both')
    if args.unmated_reads is None:
        if not all(x is not None for x in [args.mates1, args.mates2]):
            parser.error('mates are missing')
        if len(args.mates1) != len(args.mates2):
            parser.error('mismatching number of read files')

    return args

def main():
    args = parse_args()

    if not salmon_in_path():
        print('error: salmon not found in PATH; have you done module load?',
              file=sys.stderr)
        sys.exit(1)

    if args.subcommand == 'index':
        salmon_index(args.transcripts, args.index,
                     kmer_length=args.kmer_length,
                     perfect_hash=args.perfect_hash,
                     sample_interval=args.sample_interval,
                     threads=args.threads)
    elif args.subcommand == 'quant':
        salmon_quant(index=args.index,
                     read1=args.mates1 if args.mates1 is not None else args.unmated_reads,
                     read2=args.mates2, output=args.output,
                     seq_bias=args.seq_bias, gc_bias=args.gc_bias,
                     threads=args.threads, libtype=args.libtype)

if __name__ == '__main__':
    main()
