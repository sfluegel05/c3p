"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    
    A 2-oxo monocarboxylic acid typically has a carbon (second carbon) with a ketone group (=O)
    and a nearby carboxylic acid group (COOH).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")
    
    # Pattern for the 2-oxo group (C-C(=O)-R)
    oxo_keto_pattern = Chem.MolFromSmarts("[#6][#6](=O)[#6]")
    # Pattern for the carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

    oxo_matches = mol.GetSubstructMatches(oxo_keto_pattern)
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    if not oxo_matches:
        return False, "No 2-oxo group (ketone) detected"
    
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group detected"
    
    # Ensure structural integrity: proximity (two carbon atoms directly linked between groups)
    for oxo_match in oxo_matches:
        oxo_carbon_index = oxo_match[1]  # Second carbon with the keto =O group
        for carboxylic_match in carboxylic_acid_matches:
            carboxylic_index = carboxylic_match[0]  # Carbon involved in -COOH
            # Check if the oxo carbon index and carboxylic carbon index are part of the same carbon chain
            if mol.GetBondBetweenAtoms(oxo_carbon_index, carboxylic_index):
                return True, "Valid 2-oxo monocarboxylic acid structure found"
    
    return False, "Could not confirm 2-oxo monocarboxylic acid pattern in structure"