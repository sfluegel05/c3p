"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: CHEBI:26561 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is a nucleoside where the sugar component is D-ribose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the D-ribose sugar pattern with optional stereochemistry
    ribose_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](CO)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No D-ribose sugar component found"

    # Define common nucleobase patterns (including modified bases)
    nucleobase_patterns = [
        Chem.MolFromSmarts("[nH]1cnc2c1nc[nH]2"),  # Adenine
        Chem.MolFromSmarts("[nH]1c(=O)[nH]c2c1nc[nH]2"),  # Guanine
        Chem.MolFromSmarts("[nH]1ccc(N)nc1=O"),  # Cytosine
        Chem.MolFromSmarts("[nH]1ccc(=O)[nH]c1=O"),  # Uracil
        Chem.MolFromSmarts("[nH]1cnc2c1nc[nH]2[*]"),  # Modified Adenine
        Chem.MolFromSmarts("[nH]1c(=O)[nH]c2c1nc[nH]2[*]"),  # Modified Guanine
        Chem.MolFromSmarts("[nH]1ccc(N)nc1=O[*]"),  # Modified Cytosine
        Chem.MolFromSmarts("[nH]1ccc(=O)[nH]c1=O[*]"),  # Modified Uracil
    ]

    # Check if any nucleobase pattern matches
    has_nucleobase = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    if not has_nucleobase:
        return False, "No nucleobase component found"

    # Ensure the nucleobase is attached to the ribose sugar
    # The attachment point is typically the 1' position of the ribose
    attachment_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](CO)[C@H](O)[C@H]1O[*]")
    if not mol.HasSubstructMatch(attachment_pattern):
        return False, "Nucleobase not attached to the ribose sugar"

    return True, "Contains D-ribose sugar with a nucleobase attached"