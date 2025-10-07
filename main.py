from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

fasta_file='unknown_sequence.fasta'
sequences=list(SeqIO.parse(fasta_file, 'fasta'))

def gc_content(seq):
    seq=seq.upper()
    gc_count=seq.count('G')+seq.count('C')
    return round(gc_count / len(seq) * 100, 2)

gc_data=[]
for record in sequences:
    gc = gc_content(str(record.seq))
    gc_data.append({"ID": record.id, "Length": len(record.seq), "GC Content(%)": gc})

gc_df=pd.DataFrame(gc_data)
print("\nGC Content/sequence:\n", gc_df)

def motif_count(seq, motif):
    seq = seq.upper()
    motif = motif.upper()
    return seq.count(motif)

example_motif = 'ATG' #start codon by default
motif_freqs = []
for record in sequences:
    count = motif_count(str(record.seq), example_motif)
    motif_freqs.append({"ID": record.id, f"{example_motif}_count": count})

motif_df = pd.DataFrame(motif_freqs)
print(f"\nMotif ({example_motif}) counts:\n", motif_df)

#analyze codon frequency
def codon_usage(seq):
    seq = seq.upper()
    codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3])==3]
    #returns a dictionary-like counter object, keys = codons, values = counts
    return Counter(codons)

#splits
codon_counts = codon_usage(str(sequences[0].seq))
codon_df = pd.DataFrame(codon_counts.items(), columns=['Codon', 'Frequency']).sort_values(by = 'Frequency', ascending=False)

print("\nCodon frequency (first sequence):\n", codon_df.head())

plt.figure(figsize=(6, 4))
plt.bar(gc_df['ID'], gc_df['GC Content(%)'], color='teal')
plt.title('GC Content per Sequence')
plt.xlabel('Sequence ID')
plt.ylabel('GC Content (%)')
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 4))
plt.bar(codon_df["Codon"].head(10), codon_df["Frequency"].head(10), color="purple")
plt.title("Top 10 Codons (First Sequence)")
plt.xlabel("Codon")
plt.ylabel("Frequency")
plt.tight_layout()
plt.show()

