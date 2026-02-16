def brute_force(peptide, current_mass):
    global count
    if current_mass == parent_mass:
        result.append(peptide)
        count+=1
        return
    if current_mass > parent_mass:
        return
    
    for aa, mass in amino_acid_mass.items():
        brute_force(peptide + aa, current_mass + mass)

amino_acid_mass = {
    'G':57, 'A':71, 'S':87, 'P':97, 'V':99,
    'T':101, 'C':103, 'I':113, 'L':113, 'N':114,
    'D':115, 'K':128, 'Q':128, 'E':129, 'M':131,
    'H':137, 'F':147, 'R':156, 'Y':163, 'W':186
}

parent_mass = 457
global count
count = 0
result = []
brute_force("", 0)
print("Peptides with mass", parent_mass)
print(result)
print("Total peptides found:", count)