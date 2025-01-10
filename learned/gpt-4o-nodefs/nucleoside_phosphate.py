"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
from rdkit import Chem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    Nucleoside phosphates contain a nitrogenous base bonded to a sugar, with one or more phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Update phosphate group detection to include phosphate esters found commonly
    phosphate_patterns = [
        Chem.MolFromSmarts("P(=O)(O)[O-]"),  # Common phosphate group
        Chem.MolFromSmarts("P(=O)(O)O"),    # Neutral phosphate ester
    ]

    if not any(mol.HasSubstructMatch(pp) for pp in phosphate_patterns):
        return (False, "No phosphate group found")

    # Broadened SMARTS to match sugars found in various forms of configurations
    sugar_patterns = [
        Chem.MolFromSmarts("C1[C@H](O)[C@H](O)[C@@H](O)C([C@@H]1O)"),  # Ribose
        Chem.MolFromSmarts("C1[C@H](O)[C@@H](O)[C@H](O)C([C@H]1O)"),   # Deoxyribose
        Chem.MolFromSmarts("[C@@H]1O[C@H](C(O)C(O)C1)O"),             # Others in cyclic form
    ]

    if not any(mol.HasSubstructMatch(sp) for sp in sugar_patterns):
        return (False, "No recognizable sugar moiety (ribose or similar) found")

    # Updated base patterns considering flexibility
    base_patterns = [
        Chem.MolFromSmarts("c1ncnc2ncnc12"),  # Purine skeleton
        Chem.MolFromSmarts("c1n[nH]c[nH]c1"), # Pyrimidine skeleton
        Chem.MolFromSmarts("C1(C=CC=C1)[N]"), # Alterations
    ]
    
    if not any(mol.HasSubstructMatch(bp) for bp in base_patterns):
        return (False, "No recognizable nitrogenous base found")
    
    return (True, "Contains sugar, phosphate, and base consistent with a nucleoside phosphate")

# Example test case
test_smiles = "Nc1ccn([C@@H]2O[C@H](CO)[C@@H](OP(O)(O)=O)[C@H]2O)c(=O)n1"
result, reason = is_nucleoside_phosphate(test_smiles)
print(result, reason)