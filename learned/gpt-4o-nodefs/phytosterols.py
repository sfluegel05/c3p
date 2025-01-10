"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem

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

    # General steroid backbone pattern (accounting for stereochemistry and variant structures)
    steroid_backbone = Chem.MolFromSmarts("C1CCC2C1CCC3C2CCC4C3CCCC4")
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"
        
    # Check for hydroxyl group at position 3 (allowing for stereochemistry variance)
    hydroxy_position_3_pattern = Chem.MolFromSmarts("O[C@@H]1CC[C@@H](C)CC2CCCC3C2CCC4C3CCCC4")
    if not mol.HasSubstructMatch(hydroxy_position_3_pattern):
        return False, "Missing hydroxyl group at position 3"

    # Check for common unsaturations in sterols, e.g., within the B ring
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No expected double bonds found in rings"

    # Check for side chain modifications (general alkyl patterns including methyl or ethyl groups)
    side_chain_modifications = Chem.MolFromSmarts("C(C)C")
    if not mol.HasSubstructMatch(side_chain_modifications):
        return False, "No typical phytosterol side chain modifications found"

    return True, "Steroid nucleus with typical phytosterol hydroxyl, double bonds, and side chain modifications found"