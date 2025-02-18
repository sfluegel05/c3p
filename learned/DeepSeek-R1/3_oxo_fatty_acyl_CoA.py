"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: CHEBI:154527 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    Must have: 
    1. Thioester-linked CoA moiety
    2. 3-oxo group on the fatty acid (S-C(=O)-CC(=O) pattern)
    3. Adenine-containing CoA structure with phosphate groups
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for critical 3-oxo-thioester pattern: S-C(=O)-CC(=O)
    thioester_3oxo = Chem.MolFromSmarts('[S]C(=O)CC(=O)')
    if not mol.HasSubstructMatch(thioester_3oxo):
        return False, "Missing 3-oxo-thioester group (S-C(=O)-CC(=O))"
    
    # Verify CoA structure components
    # 1. Check for adenine ring pattern (n1cnc2c1ncnc2 with possible substitutions)
    adenine_pattern = Chem.MolFromSmarts('n1cnc2c([NH1,N,H0])ncnc12')
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine ring component"
    
    # 2. Check for at least two phosphate groups (OP(=O)(O)O pattern)
    phosphate_pattern = Chem.MolFromSmarts('[O]P(=O)([O])[O]')
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need ≥2"
    
    # 3. Check for pantetheine-like chain (S-C-C-N-C=O)
    pantetheine_pattern = Chem.MolFromSmarts('[S]-C-C-N-C(=O)')
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine chain (SCCNC=O)"
    
    # Optional: Basic chain length check for fatty acid part
    # Get the 3-oxo carbon chain (atoms after S-C(=O)-CC(=O))
    # This is complex; for robustness we'll skip exact length checks
    
    return True, "Contains 3-oxo-thioester, CoA structure with adenine and phosphates"