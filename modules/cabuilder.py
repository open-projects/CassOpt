import re
import os
import math

import sqlite3
import itertools

from modules import fparser

def cabuild(sqlite_file, input_file, fasta_file, pred_file, output_file, flexible_mode, n_var, hla_num):
    prediction = PredParser(pred_file)
    fasta = FastaParser(fasta_file)
    
    sb, wb = prediction.getBinders()
    hla_set = wb.keys() # getBinders() returns values for ALL hla alleles
    
    left_nodes = set()
    right_nodes = set()

    num_paths = 0

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

    #with open('pept.test', 'w') as ff:
    #    for pept in fasta.get():
    #        ff.write(pept["seq"] + "\n")
    #print("test")

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

                yield hla, binder, pept["len"], pept["seq"], pept["l_name"], pept["r_name"], pept["l_pos"], pept["r_pos"], pept["l_ins"], pept["r_ins"]
            
    cursor.executemany("INSERT INTO pept VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", pept_iter())
    conn.commit()

    cursor.execute('''CREATE TABLE stealth_junctions_pept AS
                   SELECT hla, plen, l_name, r_name, l_pos, r_pos
                   FROM pept WHERE binder = 'NB'
                   GROUP BY hla, plen, l_name, r_name, l_pos, r_pos
                   HAVING COUNT(DISTINCT l_ins || "_" || r_ins) = plen - 1 OR l_name IN ("START","start")''')
    
    cursor.execute('''SELECT COUNT(DISTINCT plen) FROM stealth_junctions_pept''')
    len_num = cursor.fetchone()
    
    cursor.execute('''CREATE TABLE stealth_junctions (
              hla TEXT,
              l_name TEXT,
              r_name TEXT,
              l_pos INT,
              r_pos INT,
              PRIMARY KEY (l_name, r_name, l_pos, r_pos, hla)
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

    cass_fasta = Cass2Fasta(input_file)
    with open(output_file, 'w') as out_file:
        #for names in permutation:
        k = 0
        m = math.factorial(len(inter_names))
        for inter in itertools.permutations(inter_names):
            k = k + 1
            names = [start_name]
            names.extend(list(inter))
            if stop_name != '':
                names.append(stop_name)
            path = Path()

            if flexible_mode:
                hla_maxset = list(hla_set)
                for i in range(len(names) - 1):
                    l_name, r_name = names[i], names[i + 1]
                    sql = "SELECT l_pos, r_pos, COUNT(DISTINCT hla) AS n FROM stealth_junctions " \
                          "WHERE l_name = ? AND r_name = ? AND hla in ({hla}) " \
                          "GROUP BY l_pos, r_pos " \
                          "ORDER BY n DESC, l_pos - r_pos LIMIT 1".format(hla=','.join(['?'] * len(hla_maxset)))
                    cursor.execute(sql, [l_name, r_name] + hla_maxset)

                    for l_pos, r_pos, n in cursor.fetchall():
                        path.append(l_name, l_pos, r_name, r_pos)
                        cursor_hla = conn.cursor()
                        sql = "SELECT DISTINCT hla FROM stealth_junctions WHERE l_name = ? AND l_pos = ? AND r_name = ? AND r_pos = ? AND hla in ({hla})".format(hla=','.join(['?'] * len(hla_maxset)))
                        cursor_hla.execute(sql, [l_name, l_pos, r_name, r_pos] + hla_maxset)
                        hla_maxset = list()
                        for hla, in cursor_hla.fetchall():
                            hla_maxset.append(hla)
                    if len(hla_maxset) == 0:
                        break

                if len(hla_maxset) > hla_num:
                    cassette_path = ''
                    for item in path.get():
                        if len(cassette_path):
                            cassette_path += '>' + str(item['l_pos']) + '|' + str(item['r_pos']) + '<' + str(item['r_name'])
                        else:
                            cassette_path = item['l_name'] + '>' + str(item['l_pos']) + '|' + str(item['r_pos']) + '<' + item['r_name']

                    num_paths += 1
                    out_file.write("{}\t{}\t{}\n".format(num_paths, ','.join(hla_maxset), cassette_path))
                    cass_fasta.write(cassette_path, output_file + '.' + str(num_paths) + '.fa')

                    if n_var > 0 and num_paths >= n_var:
                        break
            else:
                in_path = 0
                n_hla = len(hla_set)
                sql = "SELECT l_pos, r_pos, COUNT(DISTINCT hla) AS n FROM stealth_junctions " \
                      "WHERE l_name = ? AND r_name = ? " \
                      "GROUP BY l_pos, r_pos HAVING n = ?" \
                      "ORDER BY n DESC, l_pos - r_pos LIMIT 1"
                for i in range(len(names) - 1):
                    l_name, r_name =  names[i], names[i + 1]
                    cursor.execute(sql, [l_name, r_name, n_hla])
                    in_path = 0
                    for l_pos, r_pos, n in cursor.fetchall():
                        path.append(l_name, l_pos, r_name, r_pos)
                        in_path = 1
                    if not in_path:
                        break

                if in_path:
                    cassette_path = ''
                    for item in path.get():
                        if len(cassette_path):
                            cassette_path += '>' + str(item['l_pos']) + '|' + str(item['r_pos']) + '<' + str(
                                item['r_name'])
                        else:
                            cassette_path = item['l_name'] + '>' + str(item['l_pos']) + '|' + str(item['r_pos']) + '<' + \
                                            item['r_name']

                    num_paths += 1
                    out_file.write("{}\tfull_set\t{}\n".format(num_paths, cassette_path))
                    cass_fasta.write(cassette_path, output_file + '.' + str(num_paths) + '.fa')

                    if n_var > 0 and num_paths >= n_var:
                        break

            if k % 100 == 0:
                print("iteration {} of {} ({}%), found {} variants".format(k, m, int(k/m*100), num_paths), end='', flush=True)
                print('\r', end='')
    print(' ' * (len(str(k)) + len(str(m)) + len(str(num_paths)) + 40) + "\r", end='', flush=True)
    return num_paths
# end of main()

class Cass2Fasta:
    def __init__(self, input_file):
        self._fasta = fparser.FastaParser(input_file)

    def write(self, cassette, output_file):
        with open(output_file, 'w') as out_file:
            out_file.write(">" + cassette + "\n")
            for pept in cassette.split("|"):
                p_pattern = re.match("(?:(\d+)<)?([^<>]+)(?:>(\d+))?", pept)
                if p_pattern:
                    cbeg = p_pattern.group(1)
                    name = p_pattern.group(2)
                    cend = p_pattern.group(3)
                else:
                    raise Exception("Cass2Fasta: wrong cassette format")
                item = self._fasta.get(name)
                if not item:
                    raise Exception("Cass2Fasta: can't find a header in the fasta file")


                seq = item["seq"][:item["beg"]-1].lower() + item["seq"][item["beg"]-1:item["end"]].upper() + item["seq"][item["end"]:].lower()
                cbeg = int(cbeg) if cbeg else 1
                cend = int(cend) if cend else len(seq)
                cass_seq = seq[cbeg-1:cend]
                if len(cass_seq) < 1:
                    raise Exception("Cass2Fasta: a zero length peptide is found")
                out_file.write(cass_seq + "\n")
        return
# end of class Cass2Fasta

class Path:
    def __init__(self, l_name=None, l_pos=0, r_name=None, r_pos=0):
        self._chain = []
        if l_name and r_name:
            self._chain.append({'l_name': l_name, 'l_pos': l_pos, 'r_name': r_name, 'r_pos': r_pos})
    
    def append(self, l_name, l_pos, r_name, r_pos):
        if len(self._chain) == 0 or self._chain[-1]['r_name'] == l_name:
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

    def len(self):
        return len(self._chain)
# end of class PathFinder

class PredParser:
    def __init__(self, file_name):
        self._pred = file_name
        
    def getBinders(self):
        sb = {}
        wb = {}
        with open(self._pred) as pred:
            parser = re.compile(" *\d+ +(HLA-\S+) +([GPAVLIMCFYWHKRQNEDST]{8,}) +\S+ +\d+ +\d+ +\d+ +\d+ +\d+ +\S+ +\S+ +[0-9.]+ +[0-9.]+(?:.*([SW]B))?") # !!!important: 'X' MUST BE EXCLUDED, IT'S USED AS A DUMMY START
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
                    if len(h) == 7:
                        self._fasta[header] = {"l_name": h[0], "r_name": h[1], "l_pos": int(h[2]), "l_ins": int(h[3]), "r_ins": int(h[4]), "r_pos": int(h[5]), "len": int(h[6]), "seq": ""}
                    else:
                        raise ValueError('Wrong FASTA header format!')
                elif header:
                    self._fasta[header]["seq"] += re.sub("[^A-Za-z*]", "", line)
                else:
                    raise ValueError('Wrong FASTA format!')
    def get(self):
        return self._fasta.values()
# end of class FastaParser
