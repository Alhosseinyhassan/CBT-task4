              # Amino Acid Modifications #

Am = ResidueDB().getResidue("Lysine")
print(Am.getName())
print(Am.getThreeLetterCode())
print(Am.getOneLetterCode())
print(Am.getAverageWeight())
print(Am.getMonoWeight())
print(Am.getPka())
print(Am.getFormula().toString())
""""""""""""""""""""""""""""""""""""""""""""""""""""

ox = ModificationsDB().getModification("Oxidation")
print(ox.getUniModAccession())
print(ox.getUniModRecordId())
print(ox.getDiffMonoMass())
print(ox.getId())
print(ox.getFullId())
print(ox.getFullName())
print(ox.getDiffFormula())
""""""""""""""""""""""""""""""""""""""""""""""""""""

isotopes = ox.getDiffFormula().getIsotopeDistribution(CoarseIsotopePatternGenerator(5))
for iso in isotopes.getContainer():
    print(iso.getMZ(), ":", iso.getIntensity())
""""""""""""""""""""""""""""""""""""""""""""""""""""

        # Ribonucleotides #

uridine = RibonucleotideDB().getRibonucleotide(b"U")
print(uridine.getName())
print(uridine.getCode())
print(uridine.getAvgMass())
print(uridine.getMonoMass())
print(uridine.getFormula().toString())
print(uridine.isModified())
methyladenosine = RibonucleotideDB().getRibonucleotide(b"m1A")
print(methyladenosine.getName())
print(methyladenosine.isModified())
""""""""""""""""""""""""""""""""""""""""""""""""""""
# Amino Acid Sequences #

seq_3 = AASequence.fromString("DFPIANGER")
prefix = seq_3.getPrefix(4)
suffix = seq_3.getSuffix(5)
concat = seq_3 + seq_3

print("Sequence:", seq_3)
print("Prefix:", prefix)
print("Suffix:", suffix)
print("Concatenated:", concat)

mfull = seq_3.getMonoWeight()
mprecursor = seq_3.getMonoWeight(Residue.ResidueType.Full, 2)

mz = seq_3.getMonoWeight(Residue.ResidueType.Full, 2) / 2.0

print("Monoisotopic mass of peptide [M] is", mfull)
print("Monoisotopic mass of peptide precursor [M+2H]2+ is", mprecursor)
print("Monoisotopic m/z of [M+2H]2+ is", mz)
""""""""""""""""""""""""""""""""""""""""""""""""""""
print("The peptide", str(seq_3), "consists of the following amino acids:")
for aa in seq_3:
    print(aa.getName(), ":", aa.getMonoWeight())
""""""""""""""""""""""""""""""""""""""""""""""""""""
seq_2 = AASequence.fromString("C[143]PKCK(Label:13C(6)15N(2))CR")

if seq_2.hasNTerminalModification():
    print("N-Term Modification: ", seq_2.getNTerminalModification().getFullId())
if seq_2.hasCTerminalModification():
    print("C-Term Modification: ", seq_2.getCTerminalModification().getFullId())
for aa in seq_2:
    if (aa.isModified()):
        print(aa.getName(), ":", aa.getMonoWeight(), ":", aa.getModificationName())
    else:
        print(aa.getName(), ":", aa.getMonoWeight())
""""""""""""""""""""""""""""""""""""""""""""""""""""
# Molecular formula #

seq_formula = seq_3.getFormula()
print("Peptide", seq_3, "has molecular formula", seq_formula)
""""""""""""""""""""""""""""""""""""""""""""""""""""
# Isotope patterns #

coarse_isotopes = seq_formula.getIsotopeDistribution(CoarseIsotopePatternGenerator(6))
for iso in coarse_isotopes.getContainer():
    print("Isotope", iso.getMZ(), "has abundance", iso.getIntensity() * 100, "%")

fine_isotopes = seq_formula.getIsotopeDistribution(FineIsotopePatternGenerator(0.01))
for iso in fine_isotopes.getContainer():
    print("Isotope", iso.getMZ(), "has abundance", iso.getIntensity() * 100, "%")


""""""""""""""""""""""""""""""""""""""""""""""""""""
def plotIsotopeDistribution(isotope_distribution, title="Isotope distribution"):
    plt.title(title)
    distribution = {"mass": [], "abundance": []}
    for iso_2 in isotope_distribution.getContainer():
        distribution["mass"].append(iso_2.getMZ())
        distribution["abundance"].append(iso_2.getIntensity() * 100)

    bars = plt.bar(distribution["mass"], distribution["abundance"], width=0.01,
                   snap=False)

    plt.ylim([0, 110])
    plt.xticks(range(math.ceil(distribution["mass"][0]) - 2,
                     math.ceil(distribution["mass"][-1]) + 2))
    plt.xlabel("Atomic mass (u)")
    plt.ylabel("Relative abundance (%)")


plt.figure(figsize=(10, 7))
plt.subplot(1, 2, 1)
plotIsotopeDistribution(coarse_isotopes, "Isotope distribution - coarse")
plt.subplot(1, 2, 2)
plotIsotopeDistribution(fine_isotopes, "Isotope distribution - fine structure")
plt.show()
""""""""""""""""""""""""""""""""""""""""""""""""""""
# Fragment ions #

suffix = seq_3.getSuffix(3)
print("=" * 35)
print("y3 ion sequence:", suffix)
y3_formula = suffix.getFormula(Residue.ResidueType.YIon, 2)
suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0
suffix.getMonoWeight(Residue.ResidueType.XIon, 2) / 2.0
suffix.getMonoWeight(Residue.ResidueType.BIon, 2) / 2.0

print("y3 mz:", suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0)
print("y3 molecular formula:", y3_formula)
""""""""""""""""""""""""""""""""""""""""""""""""""""
# Modified Sequences #

seq_3 = AASequence.fromString("PEPTIDESEKUEM(Oxidation)CER")
print(seq_3.toUnmodifiedString())
print(seq_3.toString())
print(seq_3.toUniModString())
print(seq_3.toBracketString())
print(seq_3.toBracketString(False))
print(AASequence.fromString("DFPIAM(UniMod:35)GER"))
print(AASequence.fromString("DFPIAM[+16]GER"))
print(AASequence.fromString("DFPIAM[+15.99]GER"))
print(AASequence.fromString("DFPIAM[147]GER"))
print(AASequence.fromString("DFPIAM[147.035405]GER"))
""""""""""""""""""""""""""""""""""""""""""""""""""""
s = AASequence.fromString(".(Dimethyl)DFPIAMGER.")
print(s, s.hasNTerminalModification())
s = AASequence.fromString(".DFPIAMGER.(Label:18O(2))")
print(s, s.hasCTerminalModification())
s = AASequence.fromString(".DFPIAMGER(Phospho).")
print(s, s.hasCTerminalModification())
""""""""""""""""""""""""""""""""""""""""""""""""""""
# Proteins and FASTA files #

bsa = FASTAEntry()
bsa.sequence = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
bsa.description = "BSA Bovine Albumin (partial sequence)"
bsa.identifier = "BSA"
alb = FASTAEntry()
alb.sequence = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGE"
alb.description = "ALB Human Albumin (partial sequence)"
alb.identifier = "ALB"

entries = [bsa, alb]

f = FASTAFile()
f.store("example.fasta", entries)
entries = []
f = FASTAFile()
f.load("example.fasta", entries)
print(len(entries))
for e in entries:
    print(e.identifier, e.sequence)
