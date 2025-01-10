"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol is defined as a primary alcohol where the alcoholic hydroxy group 
    is attached to a carbon, which is itself bonded to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for an aromatic primary alcohol:
    # [$([c]C[OH])] - Aromatic atom bonded to a non-ring carbon with -OH as primary
    primary_alcohol_pattern = Chem.MolFromSmarts("[$([c]C[OH])]")

    if primary_alcohol_pattern is None:
        return False, "Invalid SMARTS pattern"

    matches = mol.GetSubstructMatches(primary_alcohol_pattern)

    if matches:
        for match in matches:
            carbon_index = match[1]  # Index of the primary carbon
            carbon_atom = mol.GetAtomWithIdx(carbon_index)
            # Ensure that the carbon is not in a larger chain (i.e., only bonded to one carbon besides the aromatic and -OH)
            bonded_carbon_atoms = [neighbor for neighbor in carbon_atom.GetNeighbors() if neighbor.GetAtomicNum() == 6]
            if len(bonded_carbon_atoms) == 1:  # Only one carbon carbon bonded, verifying primary nature
                return True, "Primary alcohol bonded to aromatic ring found"
    
    return False, "Primary alcohol not bonded to an aromatic ring"

# Example usage
# smile = "CC1=C(CO)C=CC=C1"  # 2-methylbenzyl alcohol
# print(is_aromatic_primary_alcohol(smile))