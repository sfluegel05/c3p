"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol is a glycerophospholipid with inositol esterified to the phosphate group at the sn-3 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool, str: True with reason if molecule is a glycerophosphoinositol, False otherwise with reason
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol-phosphate backbone pattern: glycerol with O-P(=O)(O)O at sn-3
    glycerol_phosphate_pattern = Chem.MolFromSmarts("C(CO)COP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol-phosphate backbone at sn-3 found"
    
    # Inositol moiety esterified to the phosphate group
    inositol_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol group found esterified to the phosphate"
    
    # Fatty acyl chains pattern: detect at least two ester linkages 
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Expected at least 2 ester linkages for fatty acids, found {len(ester_matches)}"

    return True, "Contains glycerol-phosphate backbone with inositol esterified to the phosphate group"