from Bio.Seq import Seq

text = ("GTTCCGAAAGGCTAGCGCTAGGCGCCAAGCGGCCGGTTTCCTTGGCGACGGAGAGCGCGGGAATTTTAGA"
        "TAGATTGTAATTGCGGCTGCGCGGCCGCTGCCCGTGCAGCCAGAGGATCCAGCACCTCTCTTGGGGCTTC"
        "TCCGTCCTCGGCGCTTGGAAGTACGGATCTTTTTTCTCGGAGAAAAGTTCACTGGAACTGGAAGAAATGG")

dna_seq = Seq(text)
protein_seq = str(dna_seq.translate(to_stop=False)).replace("*", "")

print("Protein Sequence:")
print(protein_seq)
print()

AA_MW = {
    'A':89,'R':174,'N':132,'D':133,'C':121,'Q':146,'E':147,'G':75,
    'H':155,'I':131,'L':131,'K':146,'M':149,'F':165,'P':115,'S':105,
    'T':119,'W':204,'Y':181,'V':117
}

k = 4

def branch_and_bound(protein, k):
    max_weight = 0
    best = ""
    
    MAX_AA = max(AA_MW.values())

    for i in range(len(protein) - k + 1):
        weight = 0

        for j in range(k):
            weight += AA_MW[protein[i + j]]

            remaining = k - j - 1

            upper_bound = weight + remaining * MAX_AA

            if upper_bound <= max_weight:
                break

        if j == k - 1 and weight > max_weight:
            max_weight = weight
            best = protein[i:i+k]

    return best, max_weight


result = branch_and_bound(protein_seq, k)
print("Branch & Bound Result:")
print("Best peptide:", result[0])
print("Weight:", result[1])
