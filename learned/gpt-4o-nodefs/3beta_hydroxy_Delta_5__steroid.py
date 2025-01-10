"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.

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

    # Core steroid structure
    steroid_core_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCC(C4)C3CC2C1")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Steroid backbone not found"

    # 3beta-hydroxy group pattern (with tolerance for stereochemistry)
    hydroxy_3beta_pattern = Chem.MolFromSmarts("[C@@H](O)[CH](C)C")
    if not mol.HasSubstructMatch(hydroxy_3beta_pattern):
        return False, "3beta-hydroxy group not found"

    # Ensuring Delta(5) double bond presence between C5 and C6
    delta5_pattern = Chem.MolFromSmarts("C=C(C)C1CCC2C(C1)CCC3C(OP)CCC4")
    c5_c6_dbond = False
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetSymbol() == 'C' and bond.GetEndAtom().GetSymbol() == 'C' and bond.GetBondTypeAsDouble() == 2.0:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            # Look for connection to typical carbon atoms of a steroid structure
            if atom1.GetDegree() >= 3 and atom2.GetDegree() >= 3:
                c5_c6_dbond = True
                break
    
    if not c5_c6_dbond:
        return False, "Delta(5) double bond not found"

    return True, "Molecule is classified as a 3beta-hydroxy-Delta(5)-steroid"