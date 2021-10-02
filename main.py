from GeneToolKit import GeneToolKit

BioTool = GeneToolKit()
BioTool.genRandomSequence(100,"DNA")

print(f"\nDNA : {BioTool.sequence}\n")

print(f"[1] + Sequence Length : {len(BioTool.sequence)} \n")
print(f"[2] + Nucleotide Frequency : {BioTool.Color(BioTool.nucleotideFrequency())} \n")
print(f"[3] + DNA/RNA Transcription : {BioTool.Color(BioTool.transcript())} \n")
print(f"[4] + DNA String + Complement + Reverse Complement :")

print(f"5' {BioTool.Color(BioTool.sequence)} 3'")
print(f"   {''.join(['|' for i in range(len(BioTool.sequence))])}   ")
print(f"3' {BioTool.Color(BioTool.complement(BioTool.sequence))} 5'")
print(f"5' {BioTool.Color(BioTool.complement(True))} 3' \n")

print(f"[5] + GC Content : {BioTool.gcContent(BioTool.seq)}%")

print(f"[6] + GC Content in subsection k=10 : {BioTool.gcContentSubsec(BioTool.sequence,k=10)} \n")

print(f"[7] + Aminoacids Sequence : {BioTool.translate_seq(BioTool.sequence)} \n")

print(f"[8] + Codon frequency (L) : {BioTool.codon_usage(BioTool.sequence,'L')} \n")

print("[9] + Reading Frames : ")
for frame in BioTool.genReadingFrames():
    print(frame)
    
print("\n[10] + All protiens in 6 open reading frames : ")
for p in BioTool.allProtiensFromOpenReadingFrames():
    print(p)