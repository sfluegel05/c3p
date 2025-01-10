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

    # Check for inositol phosphate core
    inositol_phosphate_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1OP(=O)(O)O")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "No inositol phosphate core found"

    # Check for ceramide structure (amide bond with long-chain base)
    ceramide_pattern = Chem.MolFromSmarts("N[C@@H](COP(O)(=O)O)C(=O)[C@C](O)C")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide structure detected"

    # Detect presence of additional sugar moieties, such as mannose
    additional_sugars = Chem.MolFromSmarts("OC1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O")
    if mol.HasSubstructMatch(additional_sugars):
        return True, "Contains inositol, ceramide and additional sugar components"

    # Alternatively, structure matches basic inositol phosphoceramide
    return True, "Basic inositol phosphoceramide detected with inositol and ceramide"

# Note: This function relies heavily on SMARTS patterns and may need further refinement to 
# accurately capture all variations in inositol phosphoceramides found experimentally.