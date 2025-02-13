"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid contains at least one C=C or C#C bond and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing carboxylic acid group"

    # Look for carbon-carbon double bond (C=C) or triple bond (C#C)
    unsaturation_pattern_1 = Chem.MolFromSmarts("C=C")
    unsaturation_pattern_2 = Chem.MolFromSmarts("C#C")
    
    if not mol.HasSubstructMatch(unsaturation_pattern_1) and not mol.HasSubstructMatch(unsaturation_pattern_2):
        return False, "No C=C or C#C bond found"

    return True, "Contains carboxylic acid group and C=C or C#C bond"