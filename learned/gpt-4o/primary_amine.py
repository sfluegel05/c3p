"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is characterized by the presence of an -NH2 group where
    one hydrogen atom of ammonia (NH3) is replaced by a hydrocarbyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a primary amine group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for primary amine group: Nitrogen with exactly two hydrogens and one carbon
    # Simplified SMARTS pattern to focus on the fundamental -NH2 signature
    primary_amine_pattern = Chem.MolFromSmarts("[NH2]")
    if primary_amine_pattern is not None and mol.HasSubstructMatch(primary_amine_pattern):
        # Further check to ensure the Nitrogen atom is bonded to at least one carbon
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 7:  # Nitrogen
                num_hydrogens = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
                num_carbons = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6)
                if num_hydrogens == 2 and num_carbons >= 1:
                    return True, "Primary amine group detected"
    
    return False, "No primary amine group found"