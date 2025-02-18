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
    
    # Define SMARTS patterns for key B vitamin substructures
    patterns = {
        "Thiamine (B1) core": "[NH2]C1=C(SC=C2)N=C(N)C=C1",  # Thiazole and pyrimidine parts
        "Riboflavin (B2) isoalloxazine": "O=C1N=C2C(=O)N=C(N)C2=NC1C",  # Isoalloxazine ring
        "Niacin (B3)": "C1=CC(=CN=C1)C(=O)O",  # Nicotinic acid
        "Pantothenic acid (B5)": "CC(C)(CO)C(O)C(=O)NCCC(=O)O",  # Beta-alanine + pantoic acid
        "Pyridoxine (B6)": "CC1=C(O)C(CO)=C(CO)C=N1",  # Pyridoxine structure
        "Biotin (B7)": "C1CS[C@@H]2NC(=O)N[C@H]12",  # Tetrahydrothiophene-ureido
        "Folic acid (B9) pterin": "Nc1nc2c(ncnc2n1)CNc3ccc(cc3)C(=O)",  # Pterin + PABA part
        # Cobalamin (B12) is too complex; omit for simplicity
    }
    
    for name, smarts in patterns.items():
        patt = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(patt):
            return True, f"Matches {name} pattern"
    
    # Check for pyridoxal/pyridoxamine variants (B6)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C1=NC=C(C(=O)O)C(=C1O)CO")) or \
       mol.HasSubstructMatch(Chem.MolFromSmarts("C1=NC=C(C(=O)O)C(=C1O)CN")):
        return True, "Matches pyridoxal/pyridoxamine (B6)"
    
    # Check for cobalamin (simplified)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[Co]")):  # Presence of cobalt
        return True, "Contains cobalt (B12 cobalamin)"
    
    return False, "No B vitamin substructure detected"