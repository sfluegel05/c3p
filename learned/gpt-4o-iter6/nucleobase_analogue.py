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

    # Expanded patterns for nucleobase cores
    nucleobase_core_patterns = [
        Chem.MolFromSmarts("c1[nH]cnc1"),   # Generic azole (pyrimidine-like)
        Chem.MolFromSmarts("c1nc[nH]n1"),   # Purine-like
    ]

    # Functionality-modified nucleobase patterns (e.g., extra rings, halogens, thio substitutions)
    modification_patterns = [
        Chem.MolFromSmarts("c1c[nH]c(=O)[nH]c1"),   # Uracil/thio variations
        Chem.MolFromSmarts("n1c(O)nc2[nH]ncc2n1"),   # Modifications like hydroxypurine
        Chem.MolFromSmarts("c1ncnc(O)c1"),   # Hydroxy pyrimidines
        Chem.MolFromSmarts("c1nc(=O)[nH]c1"), # Hydroxy, thio, or halo on pyrimidine
    ]
    
    # Look for nucleobase core structures
    for pattern in nucleobase_core_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Core nucleobase-like structure detected"
    
    # Check for modification patterns on a nucleobase core
    for base_pattern in nucleobase_core_patterns:
        for mod_pattern in modification_patterns:
            if mol.HasSubstructMatch(base_pattern) and mol.HasSubstructMatch(mod_pattern):
                return True, "Modified nucleobase analogue structure detected"

    # Final fallback if no clear indicator was identified
    return False, "No significant nucleobase analogue characteristics detected"