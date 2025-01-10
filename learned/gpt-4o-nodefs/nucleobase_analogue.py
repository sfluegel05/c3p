"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define updated patterns for nucleobase analogues
    nucleobase_analogue_patterns = [
        "c1ncnc2[nH]ncnc12",  # Adenine/purine core
        "c1nc[nH]c(=O)n1",  # Uracil/pyrimidine core
        "c1c[nH]c[nH]c1=O",  # Base structure of uracil and variants
        "c1nc(=O)[nH]cc1",  # Cytosine and its derivatives
        "n1c[nH]cnc1",  # Imidazole-containing ring systems
        "c1cn(cn1)[nH]c=O",  # Modifications on purine system
        "c1[nH]c(=O)nc2[nH]cnc12",  # Extended purine modifications
        "n1cnc2c1ncnc2O",  # Hydroxyl modifications on purines
        "c1c[nH]ncn1",  # Adapts to modifications around the imidazole and pyrimidine rings
    ]

    # Check for substructure matches amongst defined patterns
    for pattern in nucleobase_analogue_patterns:
        substruct = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(substruct):
            return True, f"Structure matches nucleobase analogue pattern: {pattern}"
    
    return False, "Does not match patterns or typical features of nucleobase analogues"