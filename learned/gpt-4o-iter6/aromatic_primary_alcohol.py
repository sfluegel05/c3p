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

    # Define pattern for primary alcohol attached to an aromatic ring:
    # [cR] - any aromatic atom
    # [CH2OH] - primary alcohol group
    # This pattern reflects an aromatic carbon connected to a primary alcohol group
    primary_alcohol_pattern = Chem.MolFromSmarts("[c][CX2H2][OH]")

    if primary_alcohol_pattern is None:
        return False, "Invalid SMARTS pattern"

    matches = mol.GetSubstructMatches(primary_alcohol_pattern)

    if matches:
        for match in matches:
            atom_index = match[1]  # Index of the CH2 group
            atom = mol.GetAtomWithIdx(atom_index)
            # Ensure that the CH2 group is not part of a larger chain or cycle
            # by counting hydrogens and ensuring direct bond to OH is primary
            if atom.GetTotalDegree() == 2:  # CH2 connected to only OH and one aromatic
                return True, "Primary alcohol bonded to aromatic ring found"

    return False, "Primary alcohol not bonded to an aromatic ring"

# Example usage
# smile = "CC1=C(CO)C=CC=C1"  # 2-methylbenzyl alcohol
# print(is_aromatic_primary_alcohol(smile))