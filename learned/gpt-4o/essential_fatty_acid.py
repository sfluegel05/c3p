"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are a subset of polyunsaturated fatty acids required in the diet.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a terminal carboxylic acid
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid found"
    
    # Calculate carbon count to ensure it's a reasonable length for an essential fatty acid
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 16:
        return False, f"Carbon chain is too short: only {carbon_count} carbon atoms"
    
    # Attempt to find multiple cis-configured double bonds for polyunsaturation
    # This pattern attempts to capture cis patterns along the chain
    cis_double_bond_pattern = Chem.MolFromSmarts("C/C=C\\C")
    double_bond_count = len(mol.GetSubstructMatches(cis_double_bond_pattern))
    
    # We consider a significant number of double bonds for the class
    if double_bond_count < 3:
        return False, f"Insufficient cis double bonds: only {double_bond_count} found"
    
    # Example: a fatty acid typically will not have complex ring structures or heteroatoms beyond -COOH
    # Create pattern to recognize exclusions due to complex structures, if needed
    complex_structure_pattern = Chem.MolFromSmarts("[R]")  # Generic, can refine
    if mol.HasSubstructMatch(complex_structure_pattern):
        return False, "Complex structure not typical of essential fatty acids"
    
    return True, "Matches essential fatty acid pattern with polyunsaturation and a terminal carboxylic acid"