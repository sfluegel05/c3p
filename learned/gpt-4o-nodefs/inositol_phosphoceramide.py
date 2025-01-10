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

    # Check for ceramide structure (amide bond with variable long-chain base)
    # More generalized ceramide pattern: N-[C@H](COP(O)(=O)O)[C@@H](O)C
    ceramide_pattern = Chem.MolFromSmarts("N[C@@H](CO[P](=O)(O)O)[C@@H](O)C")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide structure detected"
    
    # Validate presence of potential long hydrophobic chains typical of ceramides
    long_chain_pattern = Chem.MolFromSmarts("C(~C)(~C)(~C)")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Missing long hydrophobic chain typical of ceramide"

    # Detect presence of additional sugar moieties like mannose
    additional_sugars = Chem.MolFromSmarts("OC1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O")
    if mol.HasSubstructMatch(additional_sugars):
        return True, "Contains inositol, ceramide, and additional sugar components"

    return True, "Basic inositol phosphoceramide detected with inositol and ceramide"