"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
from rdkit import Chem

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool, str: True and reason if the molecule is an inositol phosphoceramide, False otherwise
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify inositol ring (general form)
    inositol_smarts = "C1OC(O)C(O)C(O)C(O)C1O"  # Allow variations, ignore stereochemistry
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring detected"

    # Identify phosphodiester bridge (general form)
    phospho_smarts = "OP(=O)([O-])O"  # Focus on the bridge between inositol and ceramide
    phospho_pattern = Chem.MolFromSmarts(phospho_smarts)
    if not mol.HasSubstructMatch(phospho_pattern):
        return False, "No phosphodiester bridge found"

    # Identify ceramide moiety (general form)
    ceramide_smarts = "NC(=O)C[OH]"  # More general to allow for different ceramide-like structures
    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide moiety found"

    # General checks for minimum carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, "Insufficient carbon atoms for long chain characteristic of ceramides"

    return True, "Inositol phosphoceramide structure detected"