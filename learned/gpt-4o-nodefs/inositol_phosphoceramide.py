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
        return (False, "Invalid SMILES string")

    # Check for inositol-phosphate moiety
    inositol_phosphate_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O)C1OP(=O)(O)O")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return (False, "No inositol-phosphate moiety found")

    # Check for ceramide structure (amide linkage with potential long-chain sphingoid)
    ceramide_pattern = Chem.MolFromSmarts("C(=O)N[C@@H](COP(=O)(O)O)C(O)C")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return (False, "No ceramide structure detected")
    
    # Validate presence of potential long hydrophobic chains
    long_chain_pattern = Chem.MolFromSmarts("C(CC(C)C)C")
    if mol.GetNumAtoms(mol.GetSubstructMatch(long_chain_pattern)[0]) < 20:
        return (False, "Missing or shortened long hydrophobic chain typical of ceramide")

    # Detect presence of additional sugar moieties like mannose if present
    mannose_pattern = Chem.MolFromSmarts("OC1OC(CO)C(O)C(O)C1O")
    if mol.HasSubstructMatch(mannose_pattern):
        return (True, "Contains inositol, ceramide, and mannose sugar components")

    # Confirm basic structure without additional components
    return (True, "Basic inositol phosphoceramide detected with inositol and ceramide")