"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: Methyl Sulfide
Definition: Any aliphatic sulfide (i.e. thioether where both S substituents are carbons)
in which at least one of the organyl groups attached to the sulfur is a methyl group.
To improve over the previous implementation, we impose additional conditions:
  1. The sulfur must be non–aromatic, bonded only via single bonds,
     and its both neighbors must be carbon atoms.
  2. Among its two carbon substituents, at least one must be a methyl group,
     defined as a carbon (non–aromatic) having only one heavy-atom neighbor (degree==1)
     and exactly three attached hydrogens.
  3. The other substituent (if not also a methyl group) must be aliphatic (i.e. not aromatic).
  4. To help exclude peptides or biopolymers, if the total number of amide bonds detected
     (via the SMARTS pattern "[NX3][CX3](=O)") is 2 or greater then the molecule is assumed
     to be part of a peptide or biopolymer and is not classified as a simple methyl sulfide.
"""

from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines whether the input molecule is classified as a methyl sulfide.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a methyl sulfide; False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude peptides/biopolymers; if there are 2 or more amide bonds, we assume it is not a simple small-molecule.
    amide_smarts = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if len(amide_matches) >= 2:
        return False, f"Found {len(amide_matches)} amide bonds; likely a peptide or biopolymer"
    
    # Loop over atoms: candidate sulfur atoms.
    for atom in mol.GetAtoms():
        # Check the atom is sulfur (atomic number 16)
        if atom.GetAtomicNum() != 16:
            continue
        
        # Must be non-aromatic.
        if atom.GetIsAromatic():
            continue
        
        # Must have exactly 2 bonds.
        bonds = atom.GetBonds()
        if len(bonds) != 2:
            continue
        
        # All bonds to S must be single.
        if not all(bond.GetBondType() == Chem.rdchem.BondType.SINGLE for bond in bonds):
            continue
        
        # Get the two neighbors. For an aliphatic thioether, both neighbors should be carbon.
        neighbors = atom.GetNeighbors()
        if not all(neighbor.GetAtomicNum() == 6 for neighbor in neighbors):
            continue
        
        # Identify if at least one neighbor qualifies as a methyl group.
        # A methyl group here: non-aromatic carbon with degree==1 (only bonded to S) and exactly 3 hydrogen atoms.
        methyl_found = False
        candidate_other = None  # the other substituent
        for neighbor in neighbors:
            # Check aliphatic (should already be non-aromatic) but enforce this:
            if neighbor.GetIsAromatic():
                continue
            if neighbor.GetDegree() == 1 and neighbor.GetTotalNumHs() == 3:
                methyl_found = True
            else:
                # Save the non-methyl substituent (if any)
                candidate_other = neighbor
        
        # Only if a methyl substituent is found, go on to check the other group (if it exists)
        if methyl_found:
            # If the other substituent exists and is not itself a methyl, require that it is aliphatic.
            # (It must not be aromatic.)
            if candidate_other is not None:
                if candidate_other.GetIsAromatic():
                    continue  # reject if the other substituent is aromatic.
            return True, "Found a non-aromatic thioether with a methyl substituent"
    
    return False, "No appropriate aliphatic methyl thioether substructure found"


# Example usage (you can test with the provided SMILES strings):
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
        ("S(CC[C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H]([C@H](O)C)C(O)=O)C", "Glu-Met-Thr"),
        ("CSc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O", "6-methylthioinosine"),
        ("S(CC[C@H](NC(=O)[C@@H](N)CC=1C=2C(NC1)=CC=CC2)C(=O)N[C@@H](CC(C)C)C(O)=O)C", "Trp-Met-Leu")
    ]
    
    for s, name in tests:
        result, reason = is_methyl_sulfide(s)
        print(f"SMILES: {s}\nNAME: {name}\nResult: {result}\nReason: {reason}\n")