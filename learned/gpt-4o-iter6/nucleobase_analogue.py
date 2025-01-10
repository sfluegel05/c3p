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

    # Purine and pyrimidine SMARTS patterns
    purine_pattern = Chem.MolFromSmarts("c1ncnc2n[nH]nc12")
    pyrimidine_pattern = Chem.MolFromSmarts("c1ccncn1")
    
    # Check for core nucleobase structures
    has_purine_structure = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine_structure = mol.HasSubstructMatch(pyrimidine_pattern)
    
    if has_purine_structure or has_pyrimidine_structure:
        return True, "Core nucleobase structure detected"

    # Common nucleobase analogue modifications
    modification_patterns = [
        Chem.MolFromSmarts("[OH0,S,NH1,NH2,C=O,Cl,Br,I,F]"),  # Common modifications
        Chem.MolFromSmarts("c1[nH]nnc1")  # Recognizing imidazole-like modifications
    ]

    for pattern in modification_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Nucleobase and modification patterns detected"

    return False, "No significant nucleobase analogue characteristics detected"

# Examples of how to use the function
# result, reason = is_nucleobase_analogue('SMILES_STRING_HERE')