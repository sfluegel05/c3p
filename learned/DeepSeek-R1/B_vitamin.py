"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamins (CHEBI:...)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    B vitamins include thiamine (B1), riboflavin (B2), niacin (B3), pantothenic acid (B5),
    pyridoxine (B6), biotin (B7), folic acid (B9), and cobalamin (B12).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define valid SMARTS patterns for key B vitamin substructures
    patterns = {
        "Thiamine (B1) core": "[NH2]C1=CSC=C1.[CH2]C2=NC(=NC=C2)N",  # Thiazole + pyrimidine moieties
        "Riboflavin (B2) isoalloxazine": "O=C1N=C2C(=O)N=C(N)C2=NC1C(C)(C)C",  # Isoalloxazine ring with substituents
        "Niacin (B3)": "C(=O)(O)C1=CN=CC=C1",  # Nicotinic acid (pyridine-3-carboxylic acid)
        "Pantothenic acid (B5)": "CC(C)(CO)[C@H](O)C(=O)NCCC(=O)O",  # Beta-alanine + pantoic acid
        "Pyridoxine (B6)": "C1=NC=C(C(=O)O)C(=C1O)CO",  # Pyridoxine variant
        "Biotin (B7)": "C1CS[C@H]2NC(=O)N[C@H]12",  # Biotin's tetrahydrothiophene-ureido ring
        "Folic acid (B9) pterin": "Nc1nc2c(ncnc2n1)CN",  # Pterin moiety
        "Cobalamin (B12) cobalt": "[Co]",  # Presence of cobalt atom
    }
    
    for name, smarts in patterns.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue  # Skip invalid patterns to avoid errors
        if mol.HasSubstructMatch(patt):
            return True, f"Contains {name} substructure"
    
    # Check for pyridoxal/pyridoxamine (B6) variants
    b6_patterns = [
        "[CH3]C1=C(O)C(CO)=CN=C1",  # Pyridoxine-like
        "C1=NC=C(C(=O)O)C(=C1O)CN",  # Pyridoxamine
        "C1=NC=C(C(=O)[O-])C(=C1O)CO",  # Pyridoxal phosphate variants
    ]
    for smarts in b6_patterns:
        patt = Chem.MolFromSmarts(smarts)
        if patt and mol.HasSubstructMatch(patt):
            return True, "Matches pyridoxine derivative (B6)"
    
    # Check for folate derivatives (B9)
    folate_patt = Chem.MolFromSmarts("Nc1nc2c([NH1])c(=O)[nH]c1NCC1CN(CC1)c1ccc(cc1)C(=O)")
    if folate_patt and mol.HasSubstructMatch(folate_patt):
        return True, "Contains folate structure (B9)"
    
    return False, "No characteristic B vitamin substructure detected"