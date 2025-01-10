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
    
    # Define inositol phosphate group (more specific pattern)
    inositol_phosphate_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O)C1OP(=O)(O)O")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "No inositol phosphate group found"
    
    # Define ceramide moiety pattern which includes sphingosine and variations
    ceramide_pattern = Chem.MolFromSmarts("C(=O)NC[C@@H](O)")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide moiety found"
    
    # Check for intact phosphodiester linkage between inositol and ceramide
    phosphodiester_pattern = Chem.MolFromSmarts("COP(=O)(O)O[C@H]1")
    if not mol.HasSubstructMatch(phosphodiester_pattern):
        return False, "No phosphodiester linkage found"

    # Comprehensive validation of fatty acyl chain length (R1, R2 variability)
    # Long-chain carbon pattern hinting at fatty acyl chains
    long_chain_pattern = Chem.MolFromSmarts("C[C@H](O)C(=O)")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(long_chain_matches) < 2:
        return False, "Expected more fatty acyl chain variability"

    # Confirm successful classification if all structural components are aligned
    return True, "Contains inositol phosphate group, ceramide moiety, and phosphodiester linkage"