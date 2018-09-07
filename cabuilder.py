#!/usr/bin/env python3

import re
import argparse
import os

import sqlite3
import itertools
import sys

sqlite_file = 'peptdb.sqlite'

def main():
    input_parser = argparse.ArgumentParser(description='CaBuilder: the program to build non immunogenic cassettes.')
    input_parser.add_argument('f', metavar='input_FASTA_file.fa', help='FASTA file of junction peptides; the fasta header format: >lName_rName_lPos_lIns_rIns_rPos')
    input_parser.add_argument('p', metavar='input_netMHCpan_file.txt', help='netMHCpan output file')
    
    args = input_parser.parse_args()
    fasta_file = args.f
    pred_file = args.p
    
    prediction = PredParser(pred_file)
    fasta = FastaParser(fasta_file)
    
    sb, wb = prediction.getBinders()
    hla_set = wb.keys() # getBinders() returns values for ALL hla alleles
    
    left_nodes = set()
    right_nodes = set()
    
    if os.path.isfile(sqlite_file):
        os.remove(sqlite_file)
    conn = sqlite3.connect(sqlite_file)
    cursor = conn.cursor()
    cursor.execute('''CREATE TABLE pept (
              hla TEXT,
              binder TEXT,
              plen INT,
              pept TEXT,
              l_name TEXT,
              r_name TEXT,
              l_pos INT,
              r_pos INT,
              l_ins INT,
              r_ins INT
              )''')
              
    def pept_iter():
        for hla in hla_set:
            for pept in fasta.get():
                binder = 'NB'
                if pept["seq"] in wb[hla]:
                    binder = 'WB'
                elif pept["seq"] in sb[hla]:
                    binder = 'SB'
                
                left_nodes.add(pept["l_name"])
                right_nodes.add(pept["r_name"])
                yield hla, binder, len(pept["seq"]), pept["seq"], pept["l_name"], pept["r_name"], pept["l_pos"], pept["r_pos"], pept["l_ins"], pept["r_ins"]
            
    cursor.executemany("INSERT INTO pept VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", pept_iter())
    conn.commit()
    
    start_name = ''
    for l_node in left_nodes:
        if l_node not in right_nodes:
            start_name = l_node
            break
    if start_name == '':
        exit("No Start!")
        
    stop_name = ''
    for r_node in right_nodes:
        if r_node not in left_nodes:
            stop_name = r_node
            break
    inter_names = right_nodes
    if stop_name != '':
        inter_names.discard(stop_name)
    
    permutation = [] # all possible combinations of chunks
    for inter in list(itertools.permutations(inter_names)):
        names = [start_name]
        names.extend(inter)
        if stop_name != '':
            names.append(stop_name)
        permutation.append(names)
    
    cursor.execute('''CREATE TABLE stealth_junctions_pept AS
                   SELECT hla, plen, l_name, r_name, l_pos, r_pos
                   FROM pept WHERE binder = 'NB'
                   GROUP BY hla, plen, l_name, r_name, l_pos, r_pos
                   HAVING COUNT(DISTINCT l_ins || "_" || r_ins) = plen - 1 OR l_name = "START"''')
    
    cursor.execute('''SELECT COUNT(DISTINCT plen) FROM stealth_junctions_pept''')
    len_num = cursor.fetchone()
    
    cursor.execute('''CREATE TABLE stealth_junctions (
              hla TEXT,
              l_name TEXT,
              r_name TEXT,
              l_pos INT,
              r_pos INT,
              PRIMARY KEY (hla, l_name, r_name, l_pos, r_pos)
              )''')
    cursor.execute('''INSERT INTO stealth_junctions (hla, l_name, r_name, l_pos, r_pos)
                   SELECT hla, l_name, r_name, l_pos, r_pos
                   FROM stealth_junctions_pept
                   GROUP BY hla, l_name, r_name, l_pos, r_pos
                   HAVING COUNT(DISTINCT plen) = ?''', (len_num))
    cursor.execute('''DROP TABLE pept''')
    cursor.execute('''DROP TABLE stealth_junctions_pept''')
    conn.commit()
    conn.execute("VACUUM")
    
    print("permutations: ", len(permutation))
    for hla in hla_set:
        cassettes = []
        for names in permutation:
            paths = []
            for i in range(len(names) - 1):
                l_name, r_name =  names[i], names[i + 1]
                cursor.execute('''SELECT l_pos, r_pos
                       FROM stealth_junctions
                       WHERE hla = ? AND l_name = ? AND r_name = ?''', (hla, l_name, r_name))
                path_ext = []
                for l_pos, r_pos in cursor.fetchall():
                    if i == 0:
                        path_ext.append(Path(l_name, l_pos, r_name, r_pos))
                    else:
                        for path in paths:
                            path_copy = path.copy()
                            path_copy.append(l_name, l_pos, r_name, r_pos)
                            path_ext.append(path_copy)
                print(len(path_ext), file=sys.stderr)
                paths = path_ext
                if len(paths) == 0:
                    break
            if len(paths):
                cassettes.extend(paths)
        for path in cassettes:
            casstte_path = ''
            for item in path.get():
                if len(casstte_path):
                    casstte_path += '>' + str(item['l_pos']) + '|' + str(item['r_pos']) + '<' + str(item['r_name'])
                else:
                    casstte_path = item['l_name'] + '>' + str(item['l_pos']) + '|' + str(item['r_pos']) + '<' + item['r_name']
            print(hla + "\t" + casstte_path)
# end of main()

class Path:
    def __init__(self, l_name=None, l_pos=0, r_name=None, r_pos=0):
        self._chain = []
        if l_name and r_name:
            self._chain.append({'l_name': l_name, 'l_pos': l_pos, 'r_name': r_name, 'r_pos': r_pos})
    
    def append(self, l_name, l_pos, r_name, r_pos):
        if self._chain[-1]['r_name'] == l_name:
            self._chain.append({'l_name': l_name, 'l_pos': l_pos, 'r_name': r_name, 'r_pos': r_pos})
        else:
            raise ValueError('Wrong PATH: {}>{}..{}<{}!'.format(l_name, l_pos, r_pos, r_name))
        return None
    
    def copy(self):
        new_path = Path()
        new_path.set(self.get().copy())
        return new_path
    
    def get(self):
        return self._chain
    
    def set(self, chain):
        self._chain = chain
# end of class PathFinder

class PredParser:
    def __init__(self, file_name):
        self._pred = file_name
        
    def getBinders(self):
        sb = {}
        wb = {}
        with open(self._pred) as pred:
            parser = re.compile(" *\d+ +(HLA-\S+) +([A-Z]{8,}) +\S+ +\d+ +\d+ +\d+ +\d+ +\d+ +\S+ +\S+ +[0-9.]+ +[0-9.]+(?:.*([SW]B))?")
            for line in pred:
                match = parser.match(line)
                if match:
                    hla = match.group(1)
                    peptide = match.group(2)
                    binder = match.group(3)
                    
                    if hla not in sb:
                        sb[hla] = set()
                    if hla not in wb:
                        wb[hla] = set()
                    
                    if binder == 'SB':
                        sb[hla].add(peptide)
                    elif binder == 'WB':
                        wb[hla].add(peptide)
        return (sb, wb)
# end of PredParser

class FastaParser:
    def __init__(self, file_name):
        self._fasta = {}
        with open(file_name, "r") as infile:
            header = ""
            for line in infile:
                line = line.strip()
                h_pattern = re.match(">(\S+)", line)
                if h_pattern:
                    header = h_pattern.group(1)
                    h = re.split("_", header)
                    self._fasta[header] = {"l_name": h[0], "r_name": h[1], "l_pos": int(h[2]), "l_ins": int(h[3]), "r_ins": int(h[4]), "r_pos": int(h[5]), "seq": ""}
                elif header:
                    self._fasta[header]["seq"] += re.sub("[^A-Za-z*]", "", line)
                else:
                    raise ValueError('Wrong FASTA format!')
    def get(self):
        return self._fasta.values()
# end of class FastaParser


if __name__ == "__main__":
    main()