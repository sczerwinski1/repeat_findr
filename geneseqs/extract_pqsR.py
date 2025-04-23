#!/usr/bin/env python3
from Bio import SeqIO

# 1) Read the GenBank file
gb = SeqIO.read("PAO1_PAO1.gb", "genbank")

# 2) Find the CDS with locus_tag "PA1003"
for feat in gb.features:
    if feat.type == "CDS" and feat.qualifiers.get("locus_tag", [""])[0] == "PA1003":
        cds_seq = feat.extract(gb.seq)
        gene  = feat.qualifiers.get("gene", ["pqsR"])[0]
        prod  = feat.qualifiers.get("product", [""])[0]
        header = f">{gene}|{feat.qualifiers.get('locus_tag')[0]} {prod}"
        with open("PA1003_pqsR_CDS.fasta", "w") as out:
            out.write(header + "\n")
            # wrap at 70 bp/line
            for i in range(0, len(cds_seq), 70):
                out.write(str(cds_seq[i:i+70]) + "\n")
        print("Wrote PA1003_pqsR_CDS.fasta")
        break
else:
    print("Error: CDS with locus_tag PA1003 not found in GenBank file.")
