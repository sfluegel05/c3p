"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem

def is_phosphoinositide(smiles: str):
    """
    Classifies whether a given SMILES string corresponds to a phosphoinositide.
    
    A phosphoinositide is any phosphatidylinositol with one or more phosphates on the inositol ring.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for glycerol backbone pattern with fatty acids
    glycerol_fatty_acid_pattern = Chem.MolFromSmarts("C(=O)O[CH2]C([CH2]OP(=O)(O)O)")
    if not mol.HasSubstructMatch(glycerol_fatty_acid_pattern):
        return False, "No glycerol backbone with esterified fatty acids found"
    
    # Check for inositol ring pattern
    inositol_pattern = Chem.MolFromSmarts("C1(C(C(C(C(C1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "Inositol ring not found"
    
    # Check for phosphorylated positions on inositol ring
    phospho_inositol_pattern = Chem.MolFromSmarts("C1(C(C(C(C(OP(=O)(O)O)C1O)O)O)O)O")
    phospho_inositol_matches = mol.GetSubstructMatches(phospho_inositol_pattern)
    
    if len(phospho_inositol_matches) == 0:
        return False, "No phosphorylated inositol ring found"
    
    return True, "Molecule has a phosphatidylinositol backbone with phosphorylated inositol ring"