"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol is defined as a primary alcohol where the alcoholic hydroxy group 
    is attached to a carbon which is itself bonded to an aromatic ring.

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

    # Define primary alcohol pattern (carbon attached to an OH group)
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2]O")
    
    # Find substructures matching primary alcohol
    matches = mol.GetSubstructMatches(primary_alcohol_pattern)
    
    if not matches:
        return False, "No primary alcohol found"

    # Check if attached carbon is bonded to an aromatic ring
    for match in matches:
        carbon_index = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_index)
        is_bonded_to_aromatic = any(nb.GetIsAromatic() for nb in carbon_atom.GetNeighbors())
        
        if is_bonded_to_aromatic:
            return True, "Primary alcohol bonded to aromatic ring found"
    
    return False, "Primary alcohol not bonded to an aromatic ring"

# Example usage:
# smile = "CC1=C(CO)C=CC=C1"  # 2-methylbenzyl alcohol
# print(is_aromatic_primary_alcohol(smile))