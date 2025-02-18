"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: Methyl sulfide
Definition: Any aliphatic sulfide (thioether) in which at least one of the 
organyl groups attached to the sulfur is a methyl group.
In this implementation we require that the S atom:
  1. Has exactly two neighbors and all bonds are single.
  2. Has at least one CH3 substituent (i.e. a carbon atom with only three hydrogens).
Additionally, to avoid flagging large biomolecules where a –SCH3 is only part 
of a side chain, we check for the presence of multiple amide bonds (a rough proxy
for peptides or similar molecules). If two or more amide bonds are found, the candidate
will be rejected.
"""

from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide according to the following rules:
      - The molecule must contain at least one sulfur atom (atomic number 16)
        that is involved only in single bonds and has exactly two heavy-atom neighbors.
      - Among its two neighbors, at least one must be a methyl group 
        (a carbon atom with Degree==1 and exactly three explicitly attached hydrogens).
      - To avoid flagging peptides or nucleosides (where a methylthio group is
        only a side-chain) we count amide bonds ([NX3][CX3](=O)). If 2 or more are found,
        we assume the molecule is not primarily a methyl sulfide.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a methyl sulfide, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of multiple (>=2) amide bonds.
    amide_smarts = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if len(amide_matches) >= 2:
        # Likely a peptide or a multi-residue biomolecule.
        return False, f"Found {len(amide_matches)} amide bonds; likely a peptide or biopolymer"
    
    # Loop over all atoms, searching for candidate sulfur atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue
        bonds = atom.GetBonds()
        # Must be a simple thioether: exactly 2 substituents, all bonds must be single.
        if len(bonds) != 2 or not all(bond.GetBondType() == Chem.rdchem.BondType.SINGLE for bond in bonds):
            continue
        
        neighbor_has_methyl = False
        for neighbor in atom.GetNeighbors():
            # Look for a carbon with only one heavy-atom neighbor (the S) 
            # and (if hydrogens are explicit) three hydrogens.
            if neighbor.GetAtomicNum() == 6:
                # Use explicit hydrogen count if present.
                # Note: GetTotalNumHs() returns the total number of H (implicit+explicit).
                if neighbor.GetDegree() == 1 and neighbor.GetTotalNumHs() == 3:
                    neighbor_has_methyl = True
                    break
        if neighbor_has_methyl:
            return True, "Found a sulfur atom in a thioether with a methyl substituent"
    
    return False, "No appropriate methyl thioether substructure found"

# Example usage:
if __name__ == "__main__":
    # Some test examples, including one that should be classified as a methyl sulfide:
    tests = [
        ("C(S)(=NO)CCCCCCCCSC", "9-[(methylthio)nonyl]thiohydroximic acid"),  # True positive
        ("CSCC(=O)[C@H](O)[C@H](O)COP(O)(O)=O", "1-(methylthio)ribulose 5-phosphate"),  # True positive
        ("S(CC[C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H]([C@H](O)C)C(O)=O)C", "Glu-Met-Thr"),  # False positive example (peptide)
    ]
    for s, name in tests:
        result, reason = is_methyl_sulfide(s)
        print(f"SMILES: {s}\nNAME: {name}\nResult: {result}\nReason: {reason}\n")