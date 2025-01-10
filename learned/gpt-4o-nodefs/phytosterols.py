"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols have a steroid nucleus with a hydroxyl group at the C3 position
    and often an alkyl chain modification compared to cholesterol.

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

    # General steroid backbone pattern, more flexible with stereochemistry
    steroid_backbone = Chem.MolFromSmarts("C1CC2CCC3C=C2C1CCC3")
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"

    # More generalized pattern for hydroxyl on C3
    hydroxy_position_3 = Chem.MolFromSmarts("C(C(O))CC")
    if not mol.HasSubstructMatch(hydroxy_position_3):
        return False, "Missing hydroxyl group at possible position 3"

    # Check for phytosterol side chain modifications, a general alkyl pattern
    alkyl_side_chain = Chem.MolFromSmarts("C(C)C")
    if not mol.HasSubstructMatch(alkyl_side_chain):
        return False, "No typical phytosterol side chain alkyl groups found"

    # Finally, check for the presence of double bonds within the structure
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 1:
        return False, "No expected double bonds found in structure"

    return True, "Steroid nucleus with characteristic phytosterol hydroxyl and side chain modifications found"