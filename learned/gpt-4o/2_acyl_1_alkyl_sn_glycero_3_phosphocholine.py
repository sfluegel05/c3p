"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    This class is characterized by a glycerol backbone with an alkyl ether linkage at position 1, 
    an acyl ester linkage at position 2, and a phosphocholine group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone with phosphocholine at position 3
    glycerol_phosphocholine_pattern = Chem.MolFromSmarts("[C@H](C(O[P](=O)(OCC[N+](C)(C)C)[O-])C(O)C)OC")
    if not mol.HasSubstructMatch(glycerol_phosphocholine_pattern):
        return False, "Structure doesn't match glycerol backbone with phosphocholine group"

    # Look for alkyl group at position 1
    alkyl_pattern = Chem.MolFromSmarts("O[C@H](CO)CO")
    if not mol.HasSubstructMatch(alkyl_pattern):
        return False, "No alkyl ether group found at position 1"

    # Look for acyl group at position 2
    acyl_pattern = Chem.MolFromSmarts("C(CO)(OC=O)C")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl ester group found at position 2"

    return True, "Structure matches a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, with correct glycerol backbone, alkyl, and acyl groups"