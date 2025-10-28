from alignment import align

with open("lct_exon8.txt") as f:
    lines = f.readlines()

archive = {}
next_seq = None
unknown_key = ">unknown.1_unknown_8_17"

real_names = {
    ">uc002tuu.1_hg38_8_17 1551 1 1 chr2:135808443-135809993-" : "Human",
    ">uc002tuu.1_panTro4_8_17 1551 1 1 chr2B:139763388-1397" : "Chimpanzee",
    ">uc002tuu.1_rheMac3_8_17 1551 1 1 chr13:116031545-1160" : "Rhesus Macque",
    ">uc002tuu.1_canFam3_8_17 1551 1 1 chr19:38591470-385" : "Dog",
    ">uc002tuu.1_rn5_8_17 1551 1 1 chr13:50097887-500" : "Rat",
    ">uc002tuu.1_mm10_8_17 1551 1 1 chr1:128299839-1283" : "Mouse"
}

for line in lines:
    if line[0] == '>':
        next_seq = line.strip()
    else:
        archive[next_seq] = line.strip()

for key in archive:
    if key != unknown_key:
        score, aligned_seq1, aligned_seq2 = align(archive[unknown_key], archive[key])
        print(f"Alignment score between Unknown and {real_names[key]}: {score}")
