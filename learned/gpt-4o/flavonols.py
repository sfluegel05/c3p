"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol has a core 3-hydroxyflavone structure.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated general SMARTS pattern for flavonol structure:
    # 3-hydroxyflavone structure: aromatic ring fused to a pyranone with a 3-hydroxy group
    flavonol_pattern = Chem.MolFromSmarts('Oc1c2ccccc2oc(=O)c1-c1ccc(*c1)')
    
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "No flavonol core structure found"

    # Pattern to find any hydroxyl groups, considering their common presence
    hydroxy_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)

    # Ensure at least two hydroxy groups overall, including the core one at position 3
    if len(hydroxy_matches) < 2:
        return False, f"Insufficient hydroxyl groups, found {len(hydroxy_matches)}"

    # If no other issues, confirm it's a flavonol
    return True, "Contains the core structure of a flavonol with sufficient hydroxyl groups"

# Test with some input SMILES strings as a demonstration
print(is_flavonols('O1C(C2=CC(=CC=C2)C(=O)C3=CC=CC=C3)=C(O)C(=O)C=4C1=CC(O)=CC4O'))  # Example: 2-(3-benzoylphenyl)-3,5,7-trihydroxychromen-4-one
print(is_flavonols('COc1ccc(cc1)-c1oc2cc(O)cc(O)c2c(=O)c1O'))  # Example: kaempferide