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

    # Inositol ring, allowing for substituents and flexible stereochemistry
    inositol_smarts = "C1[O][C@@H]([O])[C@H]([O])[C@@H]([O])[C@H]([O])[C@@H]1[O]"
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring detected"

    # Phosphodiester bridge pattern, more inclusive
    phospho_smarts = "O[P](=O)(O)OC" # Focus on the bridge
    phospho_pattern = Chem.MolFromSmarts(phospho_smarts)
    if not mol.HasSubstructMatch(phospho_pattern):
        return False, "No phosphodiester bridge found"

    # Ceramide moiety, flexible structure
    ceramide_smarts = "[NH][C](=O)[C@H]"
    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide moiety found"

    # Check for a long aliphatic chain linked to the ceramide
    aliphatic_chain_smarts = "C([H])[C@H](O)CCCCCCCCCC"
    long_chain_pattern = Chem.MolFromSmarts(aliphatic_chain_smarts)
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Missing long aliphatic chain typical of ceramides"

    return True, "Inositol phosphoceramide structure detected"