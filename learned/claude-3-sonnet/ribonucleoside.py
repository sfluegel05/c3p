"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: CHEBI:18254 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside consists of a nucleobase attached to a D-ribose sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - should be reasonable for a nucleoside
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:  # Most simple nucleosides are under 500 Da
        # Look for phosphate groups that would indicate a nucleotide
        phosphate_pattern = Chem.MolFromSmarts("[P](=[O])([O,OH])([O,OH])[O,OH]")
        if mol.HasSubstructMatch(phosphate_pattern):
            return False, "Contains phosphate group - likely a nucleotide"
        
        # Look for CoA-like structures
        coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
        if mol.HasSubstructMatch(coa_pattern):
            return False, "Contains CoA-like structure"

    # Look for ribose pattern with correct stereochemistry
    ribose_pattern = Chem.MolFromSmarts("[CH2]-[CH]-[CH]-[CH]-[CH]-[O]")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar found"

    # Check for hydroxyl groups on ribose
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 2:
        return False, f"Insufficient hydroxyl groups ({hydroxyl_matches})"

    # Expanded nucleobase patterns
    pyrimidine_patterns = [
        # Basic pyrimidine
        "c1[n]c[n]c1",
        # Uracil-like
        "[nX3]1c(=O)[nH]c(=O)cc1",
        # Cytosine-like
        "[nX3]1c(N)nc(=O)cc1",
        # Modified pyrimidines
        "[nX3]1c(=O)[nH]c(=S)cc1",  # thio-derivatives
        "[nX3]1c(=O)[nH]cc(F)c1",    # fluoro-derivatives
        "[nX3]1c(=O)[nH]cc(Br)c1",   # bromo-derivatives
        "[nX3]1c(=O)[nH]cc(C)c1",    # methyl-derivatives
    ]
    
    purine_patterns = [
        # Basic purine
        "c1[n]c2[n]c[n]c2[n]1",
        # Adenine-like
        "c1nc(N)[nH]c2ncnc12",
        # Guanine-like
        "c1nc(=O)[nH]c2nc(N)[nH]c12",
        # Modified purines
        "c1nc(N)[nH]c2nc(O)[nH]c12",
        "c1nc(N)nc2[nH]cnc12"
    ]

    has_base = False
    for pattern in pyrimidine_patterns + purine_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_base = True
            break
    
    if not has_base:
        return False, "No nucleobase found"

    # Check for N-glycosidic bond connecting nucleobase to ribose
    n_glycosidic_pattern = Chem.MolFromSmarts("[#7][CH]1[O][CH]([CH2][OH])[CH]([OH])[CH]1[OH]")
    if not mol.HasSubstructMatch(n_glycosidic_pattern):
        return False, "No proper N-glycosidic bond found"

    # Check ring count (should have ribose + base rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient ring count"

    # Success case
    return True, "Contains ribose sugar connected to nucleobase via N-glycosidic bond"