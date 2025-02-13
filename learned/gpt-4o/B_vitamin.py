"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamins
"""
from rdkit import Chem

def is_B_vitamin(smiles: str):
    """
    Classifies a molecule as a B vitamin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define improved and specific SMARTS patterns for B vitamins
    b_vitamins_patterns = {
        'Vitamin B1 (thiamine and derivatives)': 'Cn1cnc(C)[n+](c1C)CCOP(O)(=O)O',  # Specific to phosphorylated forms
        'Vitamin B2 (riboflavin and derivatives)': 'CNC1=NC2=C(NC1=O)c3c4c(cnc(c4)O)cnc3C2=O',  # Core isoalloxazine structure
        'Vitamin B3 (niacin and derivatives)': 'c1cc(c(nc1)C(=O)O)[Nh]',  # Add tautomeric forms
        'Vitamin B5 (pantothenic acid and derivatives)': 'NC(CCC([O-])=O)C(=O)O[C@H](CO)C(C)(C)[H]',  # Include stereochemistry
        'Vitamin B6 (pyridoxal and derivatives)': 'c1cn(C)c(CO)c(O)c1C(=O)O',  # Pyridine-based with functionalized groups
        'Vitamin B7 (biotin and derivatives)': 'OC(=O)CCC[C@H]1N2CC(C2)SC1',  # Include both ring and open chain forms
        'Vitamin B9 (folic acid and derivatives)': 'CN1C=NC2=C1c1c(c(=O)n(cn1)[nH]2)c1ccc(NC(=O)C2=NN=C3C=CC=CC3=N2)cc1',  # More detailed with folate chain
        'Vitamin B12 (cobalamin-related)': '[Co]([C@@H]1NC(CCN1C=CC2=NC[C@@H](CC2)C(C)C)C3=CNC4=C3C=CC2=C4)C',  # Focus on cobalt core structure
    }
    
    # Check for substructures
    for name, pattern in b_vitamins_patterns.items():
        substruct = Chem.MolFromSmarts(pattern)
        if not substruct:
            continue
        if mol.HasSubstructMatch(substruct):
            return True, f"Matches substructure for {name}"

    return False, "No B vitamin substructure detected"