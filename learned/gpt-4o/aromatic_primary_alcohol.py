"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has a hydroxy group attached to a carbon which is itself bonded to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a primary alcohol (-OH connected to a non-aromatic carbon)
    # which is directly bonded to an aromatic ring atom
    pattern = Chem.MolFromSmarts("[CX4H2][OH]")

    # Check for the presence of the primary alcohol pattern
    if not mol.HasSubstructMatch(pattern):
        return False, "No primary alcohol group detected"

    # SMARTS pattern to represent any aromatic atom
    aromatic_pattern = Chem.MolFromSmarts("[a]")

    # Check if any carbon bonded to -OH is also bonded to an aromatic atom
    substruct_matches = mol.GetSubstructMatches(pattern)
    for match in substruct_matches:
        carbon_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.HasSubstructMatch(aromatic_pattern):
                return True, "Contains hydroxy group attached to an aromatic ring"
        
    return False, "Alcohol group not attached to an aromatic carbon"

# The function can now be used to classify SMILES strings for aromatic primary alcohols.