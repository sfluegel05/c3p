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
        bool: True if molecule is an inositol phosphoceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify inositol ring (6-membered carbocycle with multiple hydroxyl groups)
    inositol_smarts = "C1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)O1"
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring detected"

    # Identify phosphodiester bridge
    phospho_smarts = "OP(=O)(O)OC"
    phospho_pattern = Chem.MolFromSmarts(phospho_smarts)
    if not mol.HasSubstructMatch(phospho_pattern):
        return False, "No phosphodiester bridge found"

    # Identify ceramide moiety (amide bond plus long hydrocarbon chain)
    ceramide_smarts = "N[C@@H](CO)C(=O)[C@H](O)C"
    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide moiety found"

    # Ensure sufficient length of hydrocarbon chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, "Insufficient carbon atoms for long chain characteristic of ceramides"

    return True, "Inositol phosphoceramide structure detected"