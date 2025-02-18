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
    pyridoxine (B6), biotin (B7), folate (B9), and cobalamin (B12).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define key substructures with more precise SMARTS
    patterns = {
        # Thiamine (B1) - Thiazole + aminopyrimidine connected by methylene
        "Thiamine (B1)": "[#6]C1=C(SC=[N+]1)C.[NH2]C1=NC(=NC)C=C1",
        
        # Riboflavin (B2) - Isoalloxazine ring with ribityl sidechain
        "Riboflavin (B2)": "C1(=O)NC2=C(C(=O)N1)C(=O)N=C(N2)C.C[C@H](O)[C@H](O)CO",
        
        # Niacin (B3) - Pyridine-3-carboxylic acid or amide
        "Niacin (B3)": "[#6]C(=O)O.C1=CN=CC=C1",  # Allows derivatives
        
        # Pantothenic acid (B5) - Pantoic acid + β-alanine
        "Pantothenic acid (B5)": "CC(C)(C)[C@H](O)C(=O).NCCC(=O)",  # More general
        
        # Pyridoxine (B6) - 3-hydroxy-2-methylpyridine derivatives
        "Pyridoxine (B6)": "C1=C(O)C=C(C)N=C1.[CH3].CO",  # Core + substituents
        
        # Biotin (B7) - Tetrahydrothiophene with ureido ring
        "Biotin (B7)": "C1C[S@@]C2NC(=O)N[C@H]12",
        
        # Folate (B9) - Pterin + PABA + glutamic acid
        "Folate (B9)": "Nc1nc2c(ncnc2n1)CN.c1ccc(cc1)C(=O)N",  # Pterin + PABA
        
        # Cobalamin (B12) - Corrin ring with cobalt
        "Cobalamin (B12)": "[Co]C1=C(C)C2=C3C(C=C4N5C(=CC6=[N+]C(=C(C)C=[N+]C6(C)C)[C@H](CCC(=O)N)C5(C)CC(=O)N)C(C)(C)[C@@H]4CCC(=O)N)=C(C)C1=[N+]2"
    }
    
    # Check for B12 first (complex structure)
    cobalamin_patt = Chem.MolFromSmarts(patterns["Cobalamin (B12)"])
    if cobalamin_patt and mol.HasSubstructMatch(cobalamin_patt):
        return True, "Contains cobalamin (B12) corrin structure"
    
    # Check flavin nucleotides (B2 derivatives)
    flavin_patt = Chem.MolFromSmarts("[O]P(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12")
    if flavin_patt and mol.HasSubstructMatch(flavin_patt):
        return True, "Contains flavin nucleotide (B2 derivative)"
    
    # Check other patterns
    for name, smarts in patterns.items():
        if name == "Cobalamin (B12)":  # Already checked
            continue
        patt = Chem.MolFromSmarts(smarts)
        if patt and mol.HasSubstructMatch(patt):
            return True, f"Contains {name} substructure"
    
    # Check pyridoxal phosphate variants (B6)
    plp_patt = Chem.MolFromSmarts("C1=NC=C(C(=O)[O-])C(=C1O)COP(=O)([O-])[O-]")
    if plp_patt and mol.HasSubstructMatch(plp_patt):
        return True, "Contains pyridoxal phosphate (B6)"
    
    # Check pantothenate derivatives (B5)
    pantothenate_patt = Chem.MolFromSmarts("CC(C)(C)[C@H](O)C(=O)NCCC(=O)")
    if pantothenate_patt and mol.HasSubstructMatch(pantothenate_patt):
        return True, "Contains pantothenic acid derivative (B5)"
    
    return False, "No characteristic B vitamin substructure detected"