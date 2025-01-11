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

    # Identify the inositol phosphate moiety (inositol connected to a phosphate group)
    inositol_phosphate_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O)C1OP(=O)(O)O")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return (False, "No inositol-phosphate moiety found")

    # Look for the ceramide structure (amide linkage with long chains, potential sphingoid)
    ceramide_amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(ceramide_amide_pattern):
        return (False, "No ceramide amide linkage found")

    # Check for a long chain, typical of ceramides
    # Here we just approximate a long chain with consecutive carbon atoms
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCCCCCC") # Adjust length as needed
    matches = mol.GetSubstructMatches(long_chain_pattern)
    if not matches or len(matches) < 1:
        return (False, "No long hydrophobic chain detected")

    # Detect additional mannose or sugar moieties
    mannose_pattern = Chem.MolFromSmarts("OC1OC(CO)C(O)C(O)C1O")
    if mol.HasSubstructMatch(mannose_pattern):
        return (True, "Contains inositol, ceramide, and mannose sugar components")

    # If core structures identified without additional sugars
    return (True, "Inositol phosphoceramide detected with inositol and ceramide structure")