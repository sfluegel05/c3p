"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Identifies if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is composed of a nitrogenous base, sugar, and phosphate group.

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

    # Expanded phosphate group patterns
    phosphate_patterns = [
        Chem.MolFromSmarts('[O-]P(=O)([O-])[O-]'),  # common phosphate ion form
        Chem.MolFromSmarts('OP(=O)(O)O'),           # General phosphate
        Chem.MolFromSmarts('OP(=O)([O-])O'),        # Variations with charges
        Chem.MolFromSmarts('[O-]1P([O-])([O-])O1'), # Cyclic phosphate
        Chem.MolFromSmarts('OP(=O)(O)O[P;R1]'),     # Expanded to cover connected phosphate
    ]
    
    phosphate_found = any(mol.HasSubstructMatch(pattern) for pattern in phosphate_patterns)
    if not phosphate_found:
        return False, "No phosphate group found"
    
    # Expanded sugar ring patterns (considering variations)
    ribose_patterns = [
        Chem.MolFromSmarts('OC1C(O)C(O)C(O)C1O'),   # Any-open ribose
        Chem.MolFromSmarts('C1O[C@@H]([C@H](O)[C@@H]1O)CO'),  # Specific D-ribose
        Chem.MolFromSmarts('C1O[C@H]([C@@H](O)[C@H](O)[C@H]1O)CO') # More specific stereochemistry
    ]
    deoxyribose_patterns = [
        Chem.MolFromSmarts('C1OC(CO)C(O)C1O'),      # General deoxyribose
        Chem.MolFromSmarts('C1OCC(O)C1O'),          # Specific deoxyribose pattern
    ]
    
    sugar_found = any(mol.HasSubstructMatch(pattern) for pattern in ribose_patterns + deoxyribose_patterns)
    if not sugar_found:
        return False, "No compatible sugar ring found (ribose or deoxyribose)"
    
    # Expanded nitrogenous base patterns (covering purine, pyrimidine, and variations)
    base_patterns = [
        Chem.MolFromSmarts('n1cnc2ncnc[nH]12'),       # General purine form
        Chem.MolFromSmarts('c1nc[nH]c1=O'),           # Uracil and similar
        Chem.MolFromSmarts('c1[nH]cnc1=O'),           # Thymine and similar
        Chem.MolFromSmarts('n1c[nH]cnc1'),            # Cytosine
        Chem.MolFromSmarts('n1c2c(nc[nH]c2=O)cnc1=O'), # Guanine
        Chem.MolFromSmarts('c1[cnH]c(=O)[nH]c1'),     # Variants for modified bases
    ]

    base_found = any(mol.HasSubstructMatch(pattern) for pattern in base_patterns)
    if not base_found:
        return False, "No nitrogenous base found in the structure"

    return True, "Molecule matches nucleotide structure: base, sugar, and phosphate detected"