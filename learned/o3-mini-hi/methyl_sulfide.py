"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: Methyl Sulfide
Definition: Any aliphatic sulfide (i.e. non-aromatic thioether in which S is bonded via single bonds to two carbons)
in which at least one of the organyl groups attached to the sulfur is a methyl group.
Additional requirements: if only one substituent is methyl then the other must be aliphatic.
Molecules with 2 or more amide bonds (SMARTS "[NX3][CX3](=O)") are assumed to be peptides/biopolymers and are rejected.
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines whether the input molecule is classified as a methyl sulfide.
    
    A methyl sulfide is defined as a molecule that contains at least one sulfur atom meeting the criteria:
      - The sulfur atom (atomic number 16) is non-aromatic.
      - It is bonded via two single bonds to carbon atoms only.
      - Among its two carbon substituents, at least one is a methyl group. Here a methyl group is defined as
          a non-aromatic carbon atom whose only heavy-atom neighbor is sulfur (degree==1) and that has exactly 3 hydrogens.
      - If only one substituent is a methyl group then the other substituent must be aliphatic (non–aromatic).
    Also, if two or more amide bonds are detected (via "[NX3][CX3](=O)"), the molecule is assumed to be a peptide/biopolymer.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        tuple (bool, str): (True, explanation) if the molecule is a methyl sulfide; (False, explanation) otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude peptides/biopolymers: if molecule has 2 or more amide bonds, reject.
    amide_smarts = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if len(amide_matches) >= 2:
        return False, f"Found {len(amide_matches)} amide bonds; likely a peptide or biopolymer"
    
    # Function to determine if a carbon atom qualifies as a methyl group.
    # A methyl group: non-aromatic, degree==1 (only bonded to the S) and exactly 3 attached hydrogens.
    def is_methyl(neigh):
        if neigh.GetAtomicNum() != 6:
            return False
        if neigh.GetIsAromatic():
            return False
        if neigh.GetDegree() != 1:
            return False
        if neigh.GetTotalNumHs() != 3:
            return False
        return True

    # Look for candidate sulfur atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:  # must be sulfur
            continue
        # Sulfur must be non-aromatic.
        if atom.GetIsAromatic():
            continue
        bonds = atom.GetBonds()
        if len(bonds) != 2:
            continue
        # All bonds to S must be single bonds.
        if not all(bond.GetBondType() == Chem.rdchem.BondType.SINGLE for bond in bonds):
            continue
        # Get the two neighboring atoms; they must both be carbons.
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 2 or not all(neigh.GetAtomicNum() == 6 for neigh in neighbors):
            continue

        # Evaluate each neighbor.
        # Determine if each neighbor is methyl (methyl group) and whether it is aliphatic.
        # Note: any carbon not qualifying as methyl must be aliphatic (i.e. non-aromatic) to pass.
        n1, n2 = neighbors[0], neighbors[1]
        n1_is_methyl = is_methyl(n1)
        n2_is_methyl = is_methyl(n2)
        n1_aliphatic = not n1.GetIsAromatic()
        n2_aliphatic = not n2.GetIsAromatic()
        
        # Accept the candidate if:
        #   (a) both substituents are methyl, or
        #   (b) one is methyl and the other (non-methyl) is aliphatic.
        if n1_is_methyl and n2_is_methyl:
            return True, "Found two methyl substituents on a non-aromatic thioether sulfur"
        elif n1_is_methyl and n2_aliphatic:
            return True, "Found a methyl substituent and an aliphatic substituent on the thioether sulfur"
        elif n2_is_methyl and n1_aliphatic:
            return True, "Found a methyl substituent and an aliphatic substituent on the thioether sulfur"
        # Otherwise, do not classify this sulfur as a valid methyl sulfide substructure.
        
    return False, "No appropriate aliphatic methyl thioether substructure found"


# Example usage: running some of the test cases provided.
if __name__ == "__main__":
    tests = [
        # True positives
        ("C(S)(=NO)CCCCCCCCSC", "9-[(methylthio)nonyl]thiohydroximic acid"),
        ("CSCC(=O)[C@H](O)[C@H](O)COP(O)(O)=O", "1-(methylthio)ribulose 5-phosphate"),
        ("S(CCCCCCCCN=C=S)C", "8-(methylthio)octylisothiocyanate"),
        ("CSc1ncncn1", "2-(methylthio)-1,3,5-triazine"),
        ("O=C(O)[C@@H](NO)CCCCCSC", "N-hydroxy-L-trihomomethionine"),
        ("CSc1nc(NC(C)C)nc(NC(C)C)n1", "promethin"),
        ("CSC(CC(N)C(O)=O)C(O)=O", "gamma-thiomethyl glutamate"),
        ("CSC[C@H]1O[C@H](O)[C@H](O)[C@@H]1O", "5-S-methyl-5-thio-alpha-D-ribofuranose"),
        ("CSC[C@H](NC(C)=O)C(O)=O", "N-acetyl-S-methyl-L-cysteine"),
        ("C(\\CCCCCCCCSC)=N/O", "(E)-9-(methylsulfanyl)nonanal oxime"),
        ("CSCC1=NC2=C(NC1=O)[C@@H](C)C1(C(=O)N(C)C3=C1C=CC=C3)C1=C2C(=O)CC1", "maremycin G"),
        ("CSCCCC(=O)C(O)=O", "5-methylthio-2-oxopentanoic acid"),
        ("O=C(O)[C@@H](N(O)O)CCCCSC", "N,N-dihydroxy-L-dihomomethionine"),
        ("S(CCCCCCN=C=S)C", "1-isothiocyanato-6-(methylsulfanyl)hexane"),
        ("CNC(=O)OC1=CC(C)=C(SC)C(C)=C1", "methiocarb"),
        ("C1=C(C=CC(=C1)C[C@@H](C(O)=O)OCC)OCCN2C(=CC=C2C)C=3C=CC(=CC3)SC", "saroglitazar"),
        ("C(\\CCCCCCCSC)=N/O", "(E)-8-(methylsulfanyl)octanal oxime"),
        ("[H]C(=NO)C(C)(C)SC", "2-methyl-2-(methylsulfanyl)propanal oxime"),
        ("O=C(O)[C@@H](N)CCCCCCSC", "L-tetrahomomethionine"),
        ("CCCCOC(CSC)=O", "butyl 2-(methylsulfanyl)acetate"),
        ("CSCC(O)=O", "(methylthio)acetic acid"),
        ("C(CCCCSC)=NO", "5-(methylsulfanyl)pentanal oxime"),
        # False positives (should be rejected)
        ("CSc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O", "6-methylthioinosine"),
        ("CC1=CC(=C(S1)SC)C2C(=C(N(C3=C2C(=O)CCC3)C4=CN=CC=C4)N)C#N", "2-amino-4-[5-methyl-2-(methylthio)-3-thiophenyl]-5-oxo-1-(3-pyridinyl)-4,6,7,8-tetrahydroquinoline-3-carbonitrile"),
        ("S(CCSC)C(SC)CC", "1-(Methylthio)propyl propyl disulfide"),
        ("CN1C(=O)C(=C(N=C1SC)C2=CC=C(C=C2)Cl)C#N", "4-(4-chlorophenyl)-1-methyl-2-(methylthio)-6-oxo-5-pyrimidinecarbonitrile"),
        ("S(C=1C=2C(N=C(C1)C)=CC=CC2)C", "2-Methyl-4-(methylthio)quinoline"),
        ("CCCN(CCC)S(=O)(=O)C1=CC=C(C=C1)C(=O)NC2=NN=C(O2)CSC", "4-(dipropylsulfamoyl)-N-[5-[(methylthio)methyl]-1,3,4-oxadiazol-2-yl]benzamide"),
        ("CCC1=CC=C(C=C1)C(=O)NC(CCSC)C(=O)OC", "2-[[(4-ethylphenyl)-oxomethyl]amino]-4-(methylthio)butanoic acid methyl ester"),
        ("CNC(=O)O\\N=C(\\C)C(C)SC", "butocarboxim"),
        ("S(C/C(=C/C)/C=O)C", "2-[(Methylthio)methyl]-2-butenal"),
        ("CSCC[C@H](NC([O-])=O)C([O-])=O", "N-carboxy-L-methionine(2-)"),
        ("CSC1=CC=CC=C1C(=O)C2CCCN(C2)C(=O)C3=CSN=N3", "[2-(methylthio)phenyl]-[1-[oxo(4-thiadiazolyl)methyl]-3-piperidinyl]methanone"),
        ("CSCC[C@H](NC(=O)CC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3[C@@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3CC[C@]12C)C(O)=O", "methioninoursodeoxycholic acid"),
        ("OC(CSC1=C2C(C=C(C3=CC=C(SC)C=C3)S2)=NC=N1)=O", "({6-[4-(methylsulfanyl)phenyl]thieno[3,2-d]pyrimidin-4-yl}sulfanyl)acetic acid"),
        ("CCCOP(=O)(OCCC)Oc1ccc(SC)cc1", "propaphos"),
        ("ClC(C(NC(=O)[C@H]1NC[C@@H](C1)CCC)[C@H]2O[C@H](SC)[C@H](O)[C@H]([C@H]2O)O)C", "1'-Demethylclindamycin"),
        ("CSC([C@H]([NH3+])C([O-])=O)C([O-])=O", "3-methylthioaspartate(1-)"),
        ("CSC1=CC=CC(=C1)NC(=S)N2CCCN(C2)CC3=CC=CC=C3F", "3-[(2-fluorophenyl)methyl]-N-[3-(methylthio)phenyl]-1,3-diazinane-1-carbothioamide"),
        ("C(O)C(=O)CCSC", "1-hydroxy-4-(methythio)butan-2-one"),
        ("CSCCCCCCCC([NH3+])C([O-])=O", "pentahomomethionine zwitterion"),
        ("CSCC[C@@H](N)C(O)=O", "D-methionine"),
        ("COC1=C(C=C(C=C1)C=C2C(=O)SC(=N2)SC)OC", "4-[(3,4-dimethoxyphenyl)methylidene]-2-(methylthio)-5-thiazolone"),
        ("S(CC([C@H]([C@@H](COP(=O)([O-])[O-])O)O)=O)C", "1-methylthio-D-xylulose 5-phosphate(2-)"),
        ("CSC(=S)NCc1c[nH]c2ccccc12", "brassinin"),
        ("O=C([O-])[C@@H](N(O)O)CCCCCCCCSC", "N,N-dihydroxy-L-hexahomomethioninate"),
        ("S(C=1C=C2OC(=O)C(=C2C=CC1)C=3C4=C(C=CC=C4)NC3)C", "Roseobacticide C")
    ]
    
    for s, name in tests:
        result, reason = is_methyl_sulfide(s)
        print(f"SMILES: {s}\nNAME: {name}\nResult: {result}\nReason: {reason}\n")