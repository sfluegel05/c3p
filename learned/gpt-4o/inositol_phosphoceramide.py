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
    
    # Define inositol phosphate group incorporating chirality and variability
    inositol_phosphate_pattern = Chem.MolFromSmarts("OC1[C@@H](O)C(O)C(O)C(O)[C@H]1OP(=O)(O)O")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "No inositol phosphate group found"
    
    # Define a more flexible ceramide moiety pattern including sphingoid diversity
    ceramide_pattern = Chem.MolFromSmarts("C(=O)N[C@@H]([CH2,CH][OX2])CO")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide-like moiety found"
    
    # Define phosphodiester linkage with positional variability
    phosphodiester_pattern = Chem.MolFromSmarts("COP(=O)(O)[O,*]") 
    if not mol.HasSubstructMatch(phosphodiester_pattern):
        return False, "No phosphodiester linkage found"
    
    # Long chain carbon pattern check (capturing variability in chain length)
    long_chain_pattern = Chem.MolFromSmarts("C[C@@H](O)C[C@H](=O)NC")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(long_chain_matches) < 2:
        return False, "Expected long acyl chain variations"
    
    # Checks passed successfully
    return True, "Contains inositol phosphate group, ceramide moiety, and phosphodiester linkage"