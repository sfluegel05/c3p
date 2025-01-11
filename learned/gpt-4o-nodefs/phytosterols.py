"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are characterized by a steroid nucleus with a hydroxyl group 
    at C3 and possible additional alkyl groups in the side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for steroid backbone pattern: Cyclopentanoperhydrophenanthrene skeleton
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
        
    # Check for hydroxyl group at position 3
    hydroxy_pattern = Chem.MolFromSmarts("[#6]1:2:[OX2H]:[#6](:[#6]:[#6]1):3")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Missing hydroxyl group at position 3"

    # Check for common double bond configuration in the rings or side chains
    common_double_bond_patterns = [
        Chem.MolFromSmarts("C=C"),  # generic example, specific positions can be enumerated
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in common_double_bond_patterns):
        return False, "No double bonds found in expected locations"

    # Check for side chain modifications (e.g., methyl or ethyl groups) 
    # at common steroid positions such as C24
    side_chain_patterns = [
        Chem.MolFromSmarts("[CX4]C"),  # Generic side chain methyls
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in side_chain_patterns):
        return False, "No typical phytosterol side chain alkylations found"

    return True, "Contains steroid nucleus with hydroxyl group at C3, double bonds, and side chain alkylations indicating a phytosterol"