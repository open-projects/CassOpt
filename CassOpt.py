#!/usr/bin/env python3

import os
import re
import argparse
import shutil
import string
import glob

from modules import cashuff
from modules import cabuilder


def main():
    input_parser = argparse.ArgumentParser(description='CassOpt: the program for mini-gene optimization.')
    input_parser.add_argument('-f', metavar='input_file.fa', help='FASTA file of peptides with flanks; the fasta header format: >name (beg_pept_pos..end_pept_pos)')
    input_parser.add_argument('-l', metavar='PEPTIDE_LENGTHS', nargs='+', type=int, default=[8,9,10,11], help='lengths of peptides', required=False)
    input_parser.add_argument('-m', metavar='MIN_FLANKS_LENGTH', type=int, default=10, help='min length of flanks', required=False)
    input_parser.add_argument('-a', metavar='HLA_ALLELES', nargs='+', default=['A02:01','B07:02'], help='HLA alleles', required=False)
    input_parser.add_argument('-o', metavar='/path/to/output_dir', default='output', help='path to output directory', required=False)
    input_parser.add_argument('-p', metavar='/path/to/predictor', default='netMHCpan4', help='path to a binding predictor', required=False)

    args = input_parser.parse_args()
    in_file = args.f
    pept_len = args.l
    allele_set = args.a
    min_flank = args.m
    out_dir = args.o
    predictor = args.p
    tmp_dir = out_dir + '/tmp'

    if predictor == 'netMHCpan4' and not check_tcsh():
        exit("tcsh program is not installed, please install it (netMHCpan4 uses tcsh)\n")

    if not check_predictor(predictor):
        exit(predictor + "program is not installed or misconfigured (it will be used for prediction of immunogenecity)\n")

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    print('preparing peptides for ' + predictor + ': ', end='')
    cashuff.get_pept(in_file, pept_len, min_flank, tmp_dir)
    print('Ok')

    print('binding estimation (it can take a long time):')
    alleles = 'HLA-' + ',HLA-'.join(allele_set)
    pred_output = tmp_dir + '/binding.pred'
    fasta_solid = tmp_dir + '/peptides.fasta'
    for len in pept_len:
        fasta_input = "{}/peptides.{}.fa".format(tmp_dir, len)
        command_string = "{} -l {} -a {} -f {} >> {}".format(predictor, len, alleles, fasta_input, pred_output)
        print(command_string)
        os.system(command_string)
    print('binding estimation is Ok')

    print('building cassete variants (it can take a long time):')
    sqldb = tmp_dir + '/peptdb.sqlite'
    cabuilder.cabuild(sqldb, fasta_solid, pred_output)
    print('building cassete variants is Ok')

    print('removing peptides for ' + predictor + ': ', end='')
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





