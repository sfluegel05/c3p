"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid contains a steroid backbone with a hydroxyl group at the 3beta position
    and a double bond between positions 5 and 6.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General steroid backbone pattern (tetracyclic ABCD ring system with flexibility)
    steroid_pattern = Chem.MolFromSmarts("*1CCC2C(C1)CCC3C2CCC4C3(CCC(C4)O)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # 3beta-hydroxyl group pattern
    beta_hydroxyl_pattern = Chem.MolFromSmarts("[C@H](O)C")
    if not mol.HasSubstructMatch(beta_hydroxyl_pattern):
        return False, "No 3beta-hydroxyl group found"
    
    # Delta-5 double bond pattern (double bond between the specific positions)
    delta5_double_bond_pattern = Chem.MolFromSmarts("C1=CCC2C(C1)CCC3=C2C(O)CC4C3CCCC4")
    if not mol.HasSubstructMatch(delta5_double_bond_pattern):
        return False, "No Delta-5 double bond found"
    
    return True, "Contains steroid backbone, 3beta-hydroxyl group, and Delta-5 double bond"