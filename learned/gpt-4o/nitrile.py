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
    A nitrile contains a carbon atom triple-bonded to a nitrogen atom (C#N)
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

    # Define basic nitrile pattern: carbon triple bonded to nitrogen
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")

    # Find nitrile substructures
    if not mol.HasSubstructMatch(nitrile_pattern):
        return False, "No nitrile (C#N) group found"

    # Confirm that nitrile group is not part of a more complex functional group
    matches = mol.GetSubstructMatches(nitrile_pattern)
    for match in matches:
        carbon = mol.GetAtomWithIdx(match[0])
        neighbors = carbon.GetNeighbors()
        if all(n.GetSymbol() != 'O' and n.GetSymbol() != 'N' for n in neighbors):
            return True, "Contains a nitrile (C#N) group not part of an excluded functional group"

    return False, "Nitrile (C#N) exists but as part of an excluded group (e.g., amide, cyanate)"