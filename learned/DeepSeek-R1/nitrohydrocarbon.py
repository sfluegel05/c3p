"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon where one or more hydrogens are replaced by nitro groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define nitro group pattern (Nitro group: [N+](=O)[O-])
    nitro_pattern = Chem.MolFromSmarts('[N+](=O)[O-]')
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    
    if not nitro_matches:
        return False, "No nitro groups found"
    
    # Check each nitro group is attached to a carbon atom
    for match in nitro_matches:
        nitro_n_idx = match[0]  # Index of the nitrogen in the nitro group
        nitro_atom = mol.GetAtomWithIdx(nitro_n_idx)
        # Check if the nitrogen is bonded to at least one carbon
        has_carbon_neighbor = any(neighbor.GetAtomicNum() == 6 for neighbor in nitro_atom.GetNeighbors())
        if not has_carbon_neighbor:
            return False, "Nitro group not attached to carbon"
    
    # All checks passed
    return True, "Contains nitro groups attached to hydrocarbon structure"