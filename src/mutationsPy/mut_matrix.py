#!/usr/bin/env python3

import pandas as pd
from mutationsPy.gen_context import gen_context

import argparse, sys


def sigprofiler_to_hdp(sigprofiler_path, hdp_outpath):
    sigprofiler = pd.read_csv(sigprofiler_path, sep = '\t')    
    hdp_mut_vector = gen_context(3)
    sigprofiler = sigprofiler.set_index('MutationType')
    sigprofiler.index.name = None
    hdp = sigprofiler.transpose()
    hdp = hdp[hdp_mut_vector]
    hdp.to_csv(hdp_outpath, sep = '\t', index=True, index_label=False)


def hdp_to_sigprofiler(hdp_path, sigprofiler_outpath):
    hdp = pd.read_csv(hdp_path, sep = '\t')
    sigprofiler = hdp.transpose()
    sigprofiler.sort_index(axis = 0, ascending = True, inplace = True)
    sigprofiler.rename_axis('MutationType', inplace = True)
    sigprofiler.reset_index(inplace = True)
    sigprofiler.to_csv(sigprofiler_outpath, sep = '\t', index = False)



def concat_sigprofiler_mutmat():
    pass


def get_arguments():
    parser = argparse.ArgumentParser(description='Format mutation matrix')
    subparsers = parser.add_subparsers(required=True, dest='cmd')
    
    # parser for sigprofiler_to_hdp
    parser_sigprofiler_to_hdp = subparsers.add_parser('sigprofiler_to_hdp', help='convert tab delimited mutation matrix from sigprofiler format to hdp format')
    parser_sigprofiler_to_hdp.add_argument('--sigprofiler_path', type=str, required=True, help='path to sigprofiler tab delimited file')
    parser_sigprofiler_to_hdp.add_argument('--hdp_outpath', type=str, required=True, help='path to output hdp tab delimited file')
    parser_sigprofiler_to_hdp.set_defaults(func=sigprofiler_to_hdp)
    
    # parser for hdp_to_sigprofiler
    parser_hdp_to_sigprofiler = subparsers.add_parser('hdp_to_sigprofiler', help='convert tab delimited mutation matrix from hdp format to sigprofiler format')
    parser_hdp_to_sigprofiler.add_argument('--hdp_path', type=str, required=True, help='path to hdp tab delimited file')
    parser_hdp_to_sigprofiler.add_argument('--sigprofiler_outpath', type=str, required=True, help='path to output sigprofiler tab delimited file')
    parser_hdp_to_sigprofiler.set_defaults(func=hdp_to_sigprofiler)
    
    return parser



def main():
    args = get_arguments().parse_args()
    if args.cmd == 'sigprofiler_to_hdp':
        sigprofiler_to_hdp(args.sigprofiler_path, args.hdp_outpath)
    if args.cmd == 'hdp_to_sigprofiler':
        hdp_to_sigprofiler(args.hdp_path, args.sigprofiler_outpath)
    

if __name__ == "__main__":
    sys.exit(main())