"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid contains more than one double bond in a linear aliphatic chain
    and a carboxylic acid group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group -COOH
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of double bonds in the molecule
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    
    # Polyunsaturated definition: more than one double bond
    if len(double_bond_matches) <= 1:
        return False, f"Found {len(double_bond_matches)} double bonds, need more than one for polyunsaturation"

    # Ensure the double bonds are part of non-aromatic chains
    if mol.GetNumAromaticRings() > 0:
        return False, "Contains aromatic structures, not typical for fatty acids"
    
    # Ensure linearity: more than 12 carbons in a main chain
    carbon_chain_pattern = Chem.MolFromSmarts("C" + "~C" * 11)  # At least 12 linear carbons
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Chain is too short or not linear, typical fatty acids have longer chains"

    return True, "Contains a carboxylic acid group and more than one double bond in an aliphatic chain"