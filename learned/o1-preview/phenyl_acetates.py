"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: phenyl acetate
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is an acetate ester obtained by formal condensation of the carboxy group of acetic acid with the hydroxy group of any phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for acetate ester group
    acetate_ester_pattern = Chem.MolFromSmarts('O[C](=O)C')  # Ester linkage with acetyl group
    acetate_matches = mol.GetSubstructMatches(acetate_ester_pattern)
    if not acetate_matches:
        return False, "No acetate ester group found"

    # For each acetate ester group, check if the oxygen is connected to an aromatic ring (phenol)
    for match in acetate_matches:
        oxygen_idx = match[0]  # Index of the ester oxygen atom
        oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
        # Get neighbors of the oxygen atom other than the carbonyl carbon
        neighbors = [atom for atom in oxygen_atom.GetNeighbors() if atom.GetIdx() != match[1]]
        for neighbor in neighbors:
            if neighbor.GetIsAromatic():
                return True, "Contains phenyl acetate ester linkage"
    return False, "No phenyl acetate ester linkage found"

__metadata__ = {
    'chemical_class': {
        'name': 'phenyl acetate',
        'definition': 'An acetate ester obtained by formal condensation of the carboxy group of acetic acid with the hydroxy group of any phenol.'
    }
}