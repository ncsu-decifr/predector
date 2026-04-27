#!/usr/bin/env python3
import json
import sys
import os

def read_ldjson(filepath):
    cleavages = {}
    if not os.path.exists(filepath):
        return cleavages
    with open(filepath, 'r') as f:
        for line in f:
            if not line.strip(): continue
            try:
                data = json.loads(line)
                # SignalP 6 ldjson format: match is usually "SP" or "LIPO" or "TAT" or "TATLIPO" or "OTHER"
                # The field for cleavage site stop is 'stop'
                if data.get('analysis') == 'signalp6' and data.get('match') in ['SP', 'LIPO', 'TAT', 'TATLIPO']:
                    cleavages[data['name']] = data.get('stop', 0)
            except:
                continue
    return cleavages

def read_fasta(filepath):
    sequences = []
    header = None
    seq = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences.append((header, "".join(seq)))
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if header:
            sequences.append((header, "".join(seq)))
    return sequences

def main():
    ldjson_path = sys.argv[1]
    fasta_path = sys.argv[2]
    out_all_path = sys.argv[3]      # All sequences, cleaved if SP found
    out_secreted_path = sys.argv[4] # Only secreted sequences, cleaved

    cleavages = read_ldjson(ldjson_path)
    fasta = read_fasta(fasta_path)

    with open(out_all_path, 'w') as f_all, open(out_secreted_path, 'w') as f_sec:
        for name, seq in fasta:
            if name in cleavages:
                stop = cleavages[name]
                mature_seq = seq[stop:]
                f_all.write(f">{name}\n{mature_seq}\n")
                f_sec.write(f">{name}\n{mature_seq}\n")
            else:
                f_all.write(f">{name}\n{seq}\n")

if __name__ == "__main__":
    main()
