"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is defined as any oxolane having an oxo-substituent 
    at any position on the tetrahydrofuran ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Correct pattern for a tetrahydrofuran (oxolane) ring - O1CCCC1
    oxolane_pattern = Chem.MolFromSmarts("O1CCCC1")
    oxolane_match = mol.GetSubstructMatches(oxolane_pattern)
    if not oxolane_match:
        return False, "No oxolane (tetrahydrofuran) ring found"
    
    # Search for an oxo group (C=O) directly attached to the oxolane ring
    for match in oxolane_match:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == "O":  # Looking at the oxygen in the ring
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == "C":
                        # Check for a carbonyl (C=O) at this neighbor
                        carbonyl_found = False
                        for nbor in neighbor.GetNeighbors():
                            if nbor.GetSymbol() == "O" and mol.GetBondBetweenAtoms(neighbor.GetIdx(), nbor.GetIdx()).GetBondType().name == "DOUBLE":
                                carbonyl_found = True
                        if carbonyl_found:
                            return True, "Contains oxolane (tetrahydrofuran) ring with oxo-substituent"
    
    return False, "No oxo-substituent found on the tetrahydrofuran ring"