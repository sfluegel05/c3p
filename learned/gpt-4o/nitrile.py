"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: Nitrile
"""
from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile contains a carbon atom triple-bonded to a nitrogen atom (C#N),
    which is not part of a more complex group like amides or cyanates.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define basic nitrile pattern: carbon triple-bonded to nitrogen
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    
    # Find nitrile substructures
    matches = mol.GetSubstructMatches(nitrile_pattern)
    if not matches:
        return False, "No nitrile (C#N) group found"

    # Ensure the nitrile carbon is not connected to oxygen or nitrogen, indicating complex groups
    for match in matches:
        carbon = mol.GetAtomWithIdx(match[0])
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetAtomicNum() in {7, 8} and neighbor.GetIdx() != match[1]:
                return False, "Nitrile (C#N) exists but as part of an excluded group (e.g., amide, cyanate)"

    return True, "Contains a nitrile (C#N) group not part of an excluded functional group"