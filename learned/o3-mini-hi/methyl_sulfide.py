"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: Methyl sulfide
Definition: Any aliphatic sulfide (thioether) in which at least one of the organyl groups attached to the sulfur is a methyl group.
In this improved implementation:
  1. We search for a sulfur atom (atomic number 16) that is not aromatic, has exactly two neighbors with both bonds single.
  2. At least one of the neighbors must be a methyl group (a carbon atom having degree 1 and exactly 3 hydrogens).
  3. To avoid flagging cases where the –SCH3 is just a substituent on a peptide or a very complex structure, we count amide bonds using the SMARTS [NX3][CX3](=O) as a rough proxy. Only if there are four or more amide bonds do we assume the molecule is likely a peptide or biopolymer.
"""

from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide according to the following rules:
      - The molecule must contain at least one sulfur atom (S, atomic number 16)
        that is not aromatic, is involved only in single bonds, and has exactly two neighbors.
      - Among its two neighbors, at least one must be a methyl group 
        (i.e. a carbon atom with degree==1 and exactly three attached hydrogens).
      - To reduce classifying complex molecules (e.g. peptides), we count amide bonds using the pattern [NX3][CX3](=O).
        If 4 or more amide bonds are detected, we assume the –SCH3 is only a side‐chain.
        
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a methyl sulfide, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count amide bonds to avoid flagging large peptide/polymers.
    amide_smarts = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if len(amide_matches) >= 4:
        return False, f"Found {len(amide_matches)} amide bonds; likely a peptide or biopolymer"
    
    # Loop over sulfur atoms in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue
        # Ensure the sulfur is not aromatic and has exactly 2 bonds.
        if atom.GetIsAromatic():
            continue
        bonds = atom.GetBonds()
        if len(bonds) != 2 or not all(bond.GetBondType() == Chem.rdchem.BondType.SINGLE for bond in bonds):
            continue
        
        # Look among its neighbors to see if one qualifies as a methyl group.
        for neighbor in atom.GetNeighbors():
            # Check that the neighbor is a carbon atom.
            if neighbor.GetAtomicNum() == 6:
                # We require this carbon to be aliphatic (not aromatic)
                if neighbor.GetIsAromatic():
                    continue
                # Check that the carbon has no heavy-atom neighbors other than our sulfur (i.e. degree==1)
                # and exactly 3 attached hydrogens (total number from implicit+explicit count).
                if neighbor.GetDegree() == 1 and neighbor.GetTotalNumHs() == 3:
                    return True, "Found a sulfur atom in a thioether with a methyl substituent"
    
    return False, "No appropriate aliphatic methyl thioether substructure found"

# Example usage:
if __name__ == "__main__":
    tests = [
        # True positives (should be classified as methyl sulfide)
        ("C(S)(=NO)CCCCCCCCSC", "9-[(methylthio)nonyl]thiohydroximic acid"),
        ("CSCC(=O)[C@H](O)[C@H](O)COP(O)(O)=O", "1-(methylthio)ribulose 5-phosphate"),
        ("S(CCCCCCCCN=C=S)C", "8-(methylthio)octylisothiocyanate"),
        ("CSc1ncncn1", "2-(methylthio)-1,3,5-triazine"),
        ("O=C(O)[C@@H](NO)CCCCCSC", "N-hydroxy-L-trihomomethionine"),
        # A few false positives (expected to be rejected)
        ("CSc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O", "6-methylthioinosine"),
        ("CC1=CC(=C(S1)SC)C2C(=C(N(C3=C2C(=O)CCC3)C4=CN=CC=C4)N)C#N", 
         "2-amino-4-[5-methyl-2-(methylthio)-3-thiophenyl]-5-oxo-1-(3-pyridinyl)-4,6,7,8-tetrahydroquinoline-3-carbonitrile")
    ]
    
    for s, name in tests:
        result, reason = is_methyl_sulfide(s)
        print(f"SMILES: {s}\nNAME: {name}\nResult: {result}\nReason: {reason}\n")