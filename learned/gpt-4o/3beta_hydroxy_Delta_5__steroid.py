"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid contains a 3beta-hydroxy group and a double bond between positions 5 and 6.

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
    
    # General SMARTS pattern for 3beta-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)[C;R1]")  # A 3beta-hydroxy group on a cyclic structure
    
    # General SMARTS pattern for Delta(5) double bond in steroid
    delta5_pattern = Chem.MolFromSmarts("C=C[C;R]")  # A double bond in a cyclic structure

    # Check for 3beta-hydroxy group
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3beta-hydroxy group found"
    
    # Check for Delta(5) double bond
    # We manually check by identifying positions between rings where double bonds occur
    carbon_bond_info = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                        for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2.0]
    delta_5_present = any((5 in pair and 6 in pair) or (pair == sorted(pair)[::-1] and all(5 <= idx < 7 for idx in pair)) for pair in carbon_bond_info)

    if not delta_5_present:
        return False, "No Delta(5) double bond found (C5-C6 double bond)"

    return True, "Contains 3beta-hydroxy group with Delta(5) double bond (C5-C6 double bond)"