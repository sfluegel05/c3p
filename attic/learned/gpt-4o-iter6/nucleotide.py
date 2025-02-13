"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Identifies if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate, consisting of a nitrogenous base, sugar, and phosphate group.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a nucleotide, False otherwise.
        str: Reason for classification.
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for different phosphate group occurrences
    phosphate_patterns = [
        Chem.MolFromSmarts('[OX1]P(=O)([OX1])[O-]'),  # Linear phosphates, common
        Chem.MolFromSmarts('[OX2]P(=O)([OX2])[OX2]'),  # Terminal phosphates
        Chem.MolFromSmarts('P([OX2H1])(=O)O[OX2]'),   # Cyclic phosphate
    ]
    
    phosphate_found = any(mol.HasSubstructMatch(pattern) for pattern in phosphate_patterns)
    if not phosphate_found:
        return False, "No phosphate group found"
    
    # Check for sugar rings (ribose and deoxyribose)
    ribose_patterns = [
        Chem.MolFromSmarts('C1CO[C@H](O)[C@@H]1O'),  # D-Ribose stereo and connectivity
        Chem.MolFromSmarts('C1COC(C1O)O'),           # Simplified ribose for catch-all
    ]
    deoxyribose_patterns = [
        Chem.MolFromSmarts('C1COC[C@H]1O'),          # D-Deoxyribose with connectivity
    ]

    sugar_found = any(mol.HasSubstructMatch(pattern) for pattern in ribose_patterns + deoxyribose_patterns)
    if not sugar_found:
        return False, "No compatible sugar ring found (ribose or deoxyribose)"
    
    # Check for nitrogenous bases
    base_patterns = [
        Chem.MolFromSmarts('n1cnc2c1ncnc2'),  # Purine base (common to adenine, guanine)
        Chem.MolFromSmarts('c1c[nH]c(=O)[nH]c1=O'),  # Pyrimidine base (uracil, thymine)
        Chem.MolFromSmarts('c1[nH]c(=O)[nH]c2c1ncnc2'),  # Cytosine variation
    ]

    base_found = any(mol.HasSubstructMatch(pattern) for pattern in base_patterns)
    if not base_found:
        return False, "No nitrogenous base found in the structure"

    return True, "Molecule matches nucleotide structure: base, sugar, and phosphate detected"