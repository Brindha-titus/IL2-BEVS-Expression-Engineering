#!/usr/bin/env python3 
# IL2_BEVS_Analysis.py — Screening Submission 
 
from Bio.SeqUtils.ProtParam import ProteinAnalysis 
import re 
 
raw_seq = ("MYRMQLLSCIALSLALVTNSAPTSSSTKKTQLQLEHLLLDLQMVILNGINNYKNPKLTRML" 
           "TFKFYMPKKATELKHLQCLEEELKPLEEVLNLAQSKNFHLRPRDLISNINVIVLELKG" 
           "SETTFMCEYADETATIVEFLNRWITFCQSIISTLT") 
 
anal = ProteinAnalysis(raw_seq) 
 
print("1) Physicochemical Properties") 
print(f"MW: {anal.molecular_weight():.0f} Da") 
print(f"pI: {anal.isoelectric_point():.2f}") 
print(f"GRAVY: {anal.gravy():.3f}") 
print(f"Instability Index: {anal.instability_index():.1f}") 
 
print("\n2) Signal Peptide Parsing (human consensus)") 
signal = raw_seq[:20] 
mature = raw_seq[20:] 
print(f"Signal (1–20): {signal}") 
print(f"Mature (21–): {mature}") 
 
print("\n3) Aggregation Hotspots (>=4 hydrophobic residues in a row)") 
hydro_pat = re.compile(r"[LIVFYW]{4,}") 
regions = [(m.start()+21, m.end()+20, m.group()) for m in hydro_pat.finditer(mature)] 
print("Flagged regions (mature numbering):", regions) 
 
print("\n4) Purification Strategy") 
print("Ni-NTA (IMAC) → TEV cleavage (optional) → SEC (Superdex 75)") 