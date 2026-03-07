#!/usr/bin/env python3

import sys
from collections import defaultdict

def read_fasta_lengths(fasta):
    lengths = {}
    name = None
    seq = []
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    lengths[name] = sum(len(s) for s in seq)
                name = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if name:
            lengths[name] = sum(len(s) for s in seq)
    return lengths


def parse_paf(paf_file):
    for line in open(paf_file):
        if line.startswith("#"):
            continue

        f = line.strip().split("\t")

        yield {
            "qname": f[0],
            "qlen": int(f[1]),
            "qstart": int(f[2]),
            "qend": int(f[3]),
            "strand": f[4],
            "tname": f[5],
            "tlen": int(f[6]),
            "tstart": int(f[7]),
            "tend": int(f[8]),
            "matches": int(f[9]),
            "block": int(f[10])
        }


def paf_to_chain(paf, qlens, tlens, out):
    chain_id = 1

    for aln in paf:

        score = aln["matches"]

        qname = aln["qname"]
        qsize = qlens.get(qname, aln["qlen"])
        qstart = aln["qstart"]
        qend = aln["qend"]

        tname = aln["tname"]
        tsize = tlens.get(tname, aln["tlen"])
        tstart = aln["tstart"]
        tend = aln["tend"]

        strand = aln["strand"]

        out.write(
            f"chain {score} {tname} {tsize} + {tstart} {tend} "
            f"{qname} {qsize} {strand} {qstart} {qend} {chain_id}\n"
        )

        block_size = aln["block"]

        out.write(f"{block_size}\n\n")

        chain_id += 1


def main():

    if len(sys.argv) != 5:
        print("Usage: paf2chain.py <paf> <query_fasta> <target_fasta> <out.chain>")
        sys.exit(1)

    paf_file = sys.argv[1]
    query_fa = sys.argv[2]
    target_fa = sys.argv[3]
    out_chain = sys.argv[4]

    qlens = read_fasta_lengths(query_fa)
    tlens = read_fasta_lengths(target_fa)

    paf = parse_paf(paf_file)

    with open(out_chain, "w") as out:
        paf_to_chain(paf, qlens, tlens, out)


if __name__ == "__main__":
    main()
