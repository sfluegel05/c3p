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

    # Basic purine and pyrimidine structure patterns
    purine_pattern = Chem.MolFromSmarts("c1nc2ncnc(n2c1)N")  # Generalized purine core with possible substitution
    pyrimidine_pattern = Chem.MolFromSmarts("c1ncncc1")      # Pyrimidine structure

    # Check for core nucleobase structures
    has_purine_structure = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine_structure = mol.HasSubstructMatch(pyrimidine_pattern)
    
    if has_purine_structure or has_pyrimidine_structure:
        return True, "Core nucleobase structure detected"

    # Patterns for common modifications found in nucleobase analogues
    modification_patterns = [
        Chem.MolFromSmarts("n1cnc2ncnc(n1)c2"),  # Azapurine
        Chem.MolFromSmarts("O=C1NC=2N(C=3NC=NC31)C(=O)N2"),  # Ethylene bridged
        Chem.MolFromSmarts("c1c[nH]c(=O)[nH]c1=O"),  # Thio- and halogen modifications
        Chem.MolFromSmarts("c1ncnc(c1)NC(=O)C"),  # Miscellaneous functional groups
    ]
    
    # Search patterns specific to modifications identified in nucleobase analogues
    for pattern in modification_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Modification consistent with nucleobase analogue detected"

    return False, "No significant nucleobase analogue characteristics detected"