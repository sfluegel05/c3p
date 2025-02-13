"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Look for inositol group connected to phosphate (basic pattern)
    inositol_phosphate_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O)C1OP(O)(=O)O")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "No inositol phosphate group found"

    # Look for ceramide moiety with variations
    # Sphingosine backbone: presence of hydroxy group and amide
    ceramide_pattern = Chem.MolFromSmarts("NC(=O)[C@@H](O)")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide moiety found"
    
    # Check for phosphodiester linkage
    phosphodiester_pattern = Chem.MolFromSmarts("COP(O)(=O)O")
    if not mol.HasSubstructMatch(phosphodiester_pattern):
        return False, "No phosphodiester linkage found"
      
    # Allow for variability in the fatty acyl chains (R1 and R2 substituents)
    # This is a more complex step, and in practice, may require enumeration 
    # or specific heuristics depending on expected R group structures

    # If all patterns are matched, confirm classification
    return True, "Contains inositol phosphate group, ceramide moiety, and phosphodiester linkage"