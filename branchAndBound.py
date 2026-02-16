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
    'G':57, 'A':71, 'S':87, 'P':97, 'V':99,
    'T':101, 'C':103, 'I':113, 'L':113, 'N':114,
    'D':115, 'K':128, 'Q':128, 'E':129, 'M':131,
    'H':137, 'F':147, 'R':156, 'Y':163, 'W':186
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

