#!/usr/bin/env python3

import os
import re
import argparse
import shutil
import string

from modules import cashuff
from modules import cabuilder

def main():
    input_parser = argparse.ArgumentParser(description='CassOpt: the program for optimization of mini-gene cassettes.')
    input_parser.add_argument('-f', metavar='input_file.fa', help='FASTA file of peptides with flanks; the fasta header format: >name (beg_pept_pos..end_pept_pos)', required=True)
    input_parser.add_argument('-l', metavar='PEPTIDE_LENGTHS', nargs='+', type=int, default=[8,9,10,11], help='lengths of peptides', required=False)
    input_parser.add_argument('-m', metavar='MIN_FLANKS_LENGTH', type=int, default=8, help='min length of flanks', required=False)
    input_parser.add_argument('-a', metavar='HLA_ALLELES', nargs='+', default=['A02:01','B07:02'], help='HLA alleles', required=False)
    input_parser.add_argument('-x', action='store_true', help='fleXible mode: use subset of HLA in addition to full set of HLA', required=False)
    input_parser.add_argument('-o', metavar='/path/to/output_dir', default='output', help='path to the output directory', required=False)
    input_parser.add_argument('-p', metavar='/path/to/predictor', default='netMHCpan4', help='path to the binding predictor', required=False)
    input_parser.add_argument('-n', type=int, default=0, help="print the first n variants (it doesn't guarantee the printing of optimal variants); 0 - print all found variants (default)", required=False)
    input_parser.add_argument('-k', action='store_true', help='keep temporary files intact', required=False)

    args = input_parser.parse_args()
    in_file = args.f
    pept_len = args.l
    allele_set = args.a
    flex_mode = args.x
    min_flank = args.m
    out_dir = args.o
    predictor = args.p
    n_var = args.n
    keep_tmp = args.k

    tmp_dir = out_dir + '/tmp'
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    if predictor == 'netMHCpan4' and not check_tcsh():
        exit("tcsh program is not installed, please install it (netMHCpan4 uses tcsh)\n")

    if not check_predictor(predictor):
        exit(predictor + "program is not installed or misconfigured (it will be used for prediction of immunogenecity)\n")

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    print('preparation of junction peptides for ' + predictor + ': ', end='')
    cashuff.get_pept(in_file, pept_len, min_flank, tmp_dir)
    print('Ok')

    print('finding the binders (it can take a long time):')
    alleles = 'HLA-' + ',HLA-'.join(allele_set)
    pred_output = tmp_dir + '/binding.pred'
    fasta_solid = tmp_dir + '/peptides.fasta'

    with open(fasta_solid, 'w') as fsolid_file:
        for len in pept_len:
            fasta_input = "{}/peptides.{}.fa".format(tmp_dir, len)
            with open(fasta_input, 'r') as finput_file:
                fst = finput_file.read()
                fsolid_file.write(fst)
            command_string = "{} -l {} -a {} -f {} >> {}".format(predictor, len, alleles, fasta_input, pred_output)
            print(command_string)
            os.system(command_string)
    print('strong and weak binders are found')

    print('building the cassette (it can take a long time):')
    sqldb = tmp_dir + '/peptdb.sqlite'
    cass_output = out_dir + '/cassettes.csv'
    n_path = cabuilder.cabuild(sqldb, fasta_solid, pred_output, cass_output, flex_mode, n_var)
    print('found {} cassette variants'.format(n_path))

    if not keep_tmp:
        print('clean temporary files ...')
        shutil.rmtree(tmp_dir)
    print('...done')
# end of main

def check_predictor(name):
    out_str = os.popen(name).read()
    installed = re.search('Usage.+' + name.rstrip(string.digits), out_str)
    if (installed):
        return 1
    return 0
# end of check_netMHCpan

def check_tcsh():
    out_str = os.popen('tcsh --version').read()
    installed = re.match(r'tcsh +\d', out_str)
    if (installed):
        return 1
    return 0
# end of check_tcsh

def do_netMHCpan(fasta_file):
    exit()


if __name__ == '__main__':
    main()





