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

    # Inositol ring with flexible stereochemistry
    inositol_smarts = "C1(OC(O)C(O)C(O)C(O)C1O)"  # Allow variations by ignoring stereochemistry
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring detected"

    # Phosphodiester bridge, considering more neutral states
    phospho_smarts = "OP(=O)(O)OC"  # Focus on the bridge, accepting neutral forms
    phospho_pattern = Chem.MolFromSmarts(phospho_smarts)
    if not mol.HasSubstructMatch(phospho_pattern):
        return False, "No phosphodiester bridge found"

    # Ceramide moiety, general broad pattern
    ceramide_smarts = "[NX3][CX3](=O)[C;H1,R0][O;H1,R0]"  # Allow general variations
    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide moiety found"

    # Ensure the structure has a long chain characteristic of ceramides
    aliphatic_chain_smarts = "C[C@H](O)CCCCCCCCCCCC"
    long_chain_pattern = Chem.MolFromSmarts(aliphatic_chain_smarts)
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Missing long aliphatic chain for ceramides"

    return True, "Inositol phosphoceramide structure detected"