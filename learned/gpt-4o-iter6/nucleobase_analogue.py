"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is defined as a molecule that can substitute for a normal nucleobase in nucleic acids.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Purine SMARTS pattern: two fused rings
    purine_pattern = Chem.MolFromSmarts("c1ncnc2ncnc12")
    # Pyrimidine SMARTS pattern: single six-membered ring
    pyrimidine_pattern = Chem.MolFromSmarts("c1cncnc1")
    
    # Check for core nucleobase structures
    has_purine_structure = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine_structure = mol.HasSubstructMatch(pyrimidine_pattern)
    
    if has_purine_structure or has_pyrimidine_structure:
        return True, "Core nucleobase structure detected"

    # Common nucleobase analogue modifications: extra rings, unusual groups, etc.
    modification_patterns = [
        Chem.MolFromSmarts("Nc1ncnc2ncnc12"),  # 8-aza modifications
        Chem.MolFromSmarts("O=C1NC(=O)C=C1"),  # 5-fluoro analogues
        Chem.MolFromSmarts("O=C1C=CC(=O)NC1"),  # Thioanalogues
        Chem.MolFromSmarts("c1n[nH]c(=O)c2c1ncn2"),  # Imidazole incorporation
    ]
    
    for pattern in modification_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Modification consistent with nucleobase analogue detected"

    return False, "No significant nucleobase analogue characteristics detected"

# Test the function on some sample data
# result, reason = is_nucleobase_analogue('SMILES_STRING_HERE')