"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime contains an aliphatic chain and an aldehyde group
    converted to an oxime (R-C=NOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for the classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the oxime group pattern "-C=NO"
    oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2][OH]")
    oxime_matches = mol.GetSubstructMatches(oxime_pattern)
    if not oxime_matches:
        return False, "No oxime group found"

    # Confirm the presence of an aliphatic chain
    aliphatic_chain_pattern = Chem.MolFromSmarts("[CX4H3,CX4H2!R]~[CX4H2,CX4H,!R]~[CX4H2,CX3H,!R]")
    aliphatic_matches = mol.GetSubstructMatches(aliphatic_chain_pattern)
    
    if not aliphatic_matches:
        return False, "No aliphatic chain found"

    for match in oxime_matches:
        carbon_atom = match[0]
        # Check if the carbon is part of an aliphatic (non-aromatic) environment
        atom = mol.GetAtomWithIdx(carbon_atom)
        if atom.GetIsAromatic() or atom.IsInRing():
            return False, "Oxime carbon is part of an aromatic ring or cyclic structure"

        # Check all atoms in the substructure for non-aromatic character
        aliphatic = True
        for atom_idx in match:
            if mol.GetAtomWithIdx(atom_idx).GetIsAromatic() or mol.GetAtomWithIdx(atom_idx).IsInRing():
                aliphatic = False
                break
        if aliphatic:
            return True, "Contains an aliphatic aldoxime group"

    return False, "Oxime group is not part of a purely aliphatic chain"