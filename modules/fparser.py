import re

class FastaParser:
    def __init__(self, file_name):
        self._fasta = {}
        with open(file_name, "r") as infile:
            header = ""
            for line in infile:
                line = line.strip()
                h_pattern = re.match(">(\S+)\D+(\d+)\.\.(\d+)", line)
                if h_pattern:
                    header = h_pattern.group(1)
                    self._fasta[header] = {"name": header, "beg": int(h_pattern.group(2)), "end": int(h_pattern.group(3)), "seq": ""}
                elif header:
                    self._fasta[header]["seq"] += re.sub("[^A-Za-z*]", "", line)
                else:
                    raise ValueError('Wrong FASTA format!')
    def get(self, header=None):
        if header:
            if header in self._fasta:
                return self._fasta[header]
        else:
            return self._fasta.values()
        return None
# end of class FastaParser