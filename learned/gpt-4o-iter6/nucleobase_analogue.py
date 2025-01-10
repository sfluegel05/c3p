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

    # Refined patterns for nucleobase cores and modifications
    nucleobase_core_patterns = [
        Chem.MolFromSmarts("c1ncnc2ncnn12"),  # Purine core
        Chem.MolFromSmarts("c1[nH]cnc2[nH]ncc12"),  # Modified purine core
        Chem.MolFromSmarts("c1ncnc1=O"),  # Pyrimidine core with mod
        Chem.MolFromSmarts("c1nc[nH]c(=O)n1"),  # Uracil and similar
    ]

    # Broadening modification to include known analogues alterations
    modification_patterns = [
        Chem.MolFromSmarts("n1c(O)nc2ncnn12"),  # Hydroxypurine
        Chem.MolFromSmarts("c1nc(=O)[nH]c(N)c1"),  # N-alkylated uracil
        Chem.MolFromSmarts("c1nc(=S)[nH]c1"),  # Thiouracil
        Chem.MolFromSmarts("c1ncnc2n[nH]nc12"),  # Azapurine variants
    ]
    
    # Look for nucleobase core structures
    for pattern in nucleobase_core_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Core nucleobase-like structure detected"
    
    # Check for modification patterns on a nucleobase core
    for mod_pattern in modification_patterns:
        if mol.HasSubstructMatch(mod_pattern):
            return True, "Modified nucleobase analogue structure detected"

    # If neither core nor modification patterns match
    return False, "No significant nucleobase analogue characteristics detected"