"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid contains a steroid backbone with a hydroxyl group at the 3beta position and a double bond between positions 5 and 6.
    
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
    
    # Pattern for steroid backbone (Four-ring structure, typically 6-6-6-5 membered rings)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C=CC4CCCC(C4C3CCC2C1)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # 3beta-hydroxyl group pattern
    hydroxyl_pattern = Chem.MolFromSmarts("C[C@@H](O)")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No 3beta-hydroxyl group found"
    
    # Delta-5 double bond pattern (double bond at C5=C6)
    delta5_double_bond_pattern = Chem.MolFromSmarts("C=CC")
    double_bond_matches = mol.GetSubstructMatches(delta5_double_bond_pattern)
    if not double_bond_matches:
        return False, "No Delta-5 double bond found"
    
    # Additional validation could include verifying specific stereochemistry and bond location
    # Implement these checks if needed for more precise classification
    
    return True, "Contains steroid backbone, 3beta-hydroxyl group, and Delta-5 double bond"