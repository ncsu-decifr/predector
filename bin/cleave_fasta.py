#!/usr/bin/env python3
import json
import sys
import os
import re

def read_ldjson(filepath):
    cleavages = {}
    if not os.path.exists(filepath):
        print(f"Warning: {filepath} not found.")
        return cleavages
    
    # Regex to extract the first number from "CS pos: 20-21. Pr: 0.9763"
    cs_regex = re.compile(r"CS pos: (\d+)-")

    with open(filepath, 'r') as f:
        for line in f:
            if not line.strip(): continue
            try:
                data = json.loads(line)
                if data.get('analysis') == 'signalp6':
                    inner_data = data.get('data', {})
                    prediction = inner_data.get('prediction', '')
                    name = inner_data.get('name', data.get('name')) # Fallback to top level name
                    
                    if prediction in ['SP', 'LIPO', 'TAT', 'TATLIPO']:
                        cs_pos = inner_data.get('cs_pos', '')
                        stop = 0
                        if cs_pos:
                            match = cs_regex.search(cs_pos)
                            if match:
                                stop = int(match.group(1))
                        
                        if name:
                            cleavages[name] = stop
            except Exception as e:
                print(f"Error parsing line: {e}")
                continue
    return cleavages

def read_fasta(filepath):
    sequences = []
    header = None
    seq = []
    if not os.path.exists(filepath):
        print(f"Warning: {filepath} not found.")
        return sequences
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
    if len(sys.argv) < 5:
        print("Usage: cleave_fasta.py <ldjson> <fasta> <out_all> <out_secreted>")
        sys.exit(1)

    ldjson_path = sys.argv[1]
    fasta_path = sys.argv[2]
    out_all_path = sys.argv[3]      # All sequences, cleaved if SP found
    out_secreted_path = sys.argv[4] # Only secreted sequences, cleaved

    cleavages = read_ldjson(ldjson_path)
    fasta = read_fasta(fasta_path)

    print(f"Loaded {len(cleavages)} SignalP 6 matches from {ldjson_path}")
    print(f"Loaded {len(fasta)} sequences from {fasta_path}")

    secreted_count = 0
    with open(out_all_path, 'w') as f_all, open(out_secreted_path, 'w') as f_sec:
        for name, seq in fasta:
            if name in cleavages:
                secreted_count += 1
                stop = cleavages[name]
                mature_seq = seq[stop:]
                f_all.write(f">{name}\n{mature_seq}\n")
                f_sec.write(f">{name}\n{mature_seq}\n")
            else:
                f_all.write(f">{name}\n{seq}\n")
    
    print(f"Wrote {len(fasta)} sequences to {out_all_path} (cleaved if SP found)")
    print(f"Wrote {secreted_count} mature secreted sequences to {out_secreted_path}")

if __name__ == "__main__":
    main()
