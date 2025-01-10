"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
from rdkit import Chem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate has a nucleobase, a sugar moiety, and one or more phosphate groups attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Extended nucleobase patterns (covering purine and pyrimidine with common modifications)
    nucleobase_patterns = [
        Chem.MolFromSmarts("c1[nH]cnc2c1ncnc2"),   # Purines (adenine-like)
        Chem.MolFromSmarts("n1cnc2c1ncnc2"),      # Purines (guaninie-like)
        Chem.MolFromSmarts("c1ccncn1"),           # Pyrimidine-like
        Chem.MolFromSmarts("c1nc[nH]c1=O"),       # Pyrimidone-like
        Chem.MolFromSmarts("n1c(=O)c[nH]cn1"),    # Cytosine-like
        Chem.MolFromSmarts("c1ccc(nc1)=O")        # Modified pyrimidines
    ]
    
    nucleobase_found = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    if not nucleobase_found:
        return False, "No nucleobase found"
    
    # Refined sugar moiety detection
    sugar_patterns = [
        Chem.MolFromSmarts("C1OC(CO)C(O)C1"),  # Ribose / Deoxyribose and basic modifications considered
        Chem.MolFromSmarts("C1OCC(O)C(O)C1"),  # Alternate sugar backbone for modifications
        Chem.MolFromSmarts("C1OCC(C)C1")       # Tetrahydrofuran-like (ribose precursor)
    ]
    
    sugar_found = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not sugar_found:
        return False, "No sugar moiety found"
    
    # Flexible phosphate group detection, allowing for mono- and polyphosphate
    phosphate_patterns = [
        Chem.MolFromSmarts("P(=O)(O)O"),  # Phosphate
        Chem.MolFromSmarts("P(=O)(O)(O)O"),  # Diphosphate link possibility, if misunderstood as poly
        Chem.MolFromSmarts("OP(O)(=O)OP(=O)(O)O") # Triphosphate (ATP/GTP like)
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in phosphate_patterns):
        return False, "No phosphate group found"

    return True, "Contains nucleobase, sugar moiety, and phosphate group"