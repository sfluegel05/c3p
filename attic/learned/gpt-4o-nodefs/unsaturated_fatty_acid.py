"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    Unsaturated fatty acids have a carboxylic acid group, at least one double bond, and a long hydrocarbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a carboxylic acid group (COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for the presence of at least one double bond (C=C)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No double bonds found, indicating a saturated fatty acid"

    # Count the number of carbons to ensure a long hydrocarbon chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, f"Too few carbon atoms ({carbon_count}), typically a shorter chain"
    
    if carbon_count > 24:
        return False, f"Too many carbon atoms ({carbon_count}), typically longer than expected for common fatty acids"

    return True, "Contains carboxylic acid group, double bonds, and suitable carbon chain length characteristic of unsaturated fatty acids"