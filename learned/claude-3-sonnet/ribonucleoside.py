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

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:  # Increased threshold for complex modifications
        return False, "Molecular weight too high"

    # Look for phosphate groups that would indicate a nucleotide
    phosphate_pattern = Chem.MolFromSmarts("[P](=[O])([O,OH])([O,OH])[O,OH]")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Contains phosphate group - likely a nucleotide"

    # Check for ribose with correct stereochemistry
    # Beta-D-ribofuranose pattern with explicit stereochemistry
    ribose_pattern = Chem.MolFromSmarts("[CH2][C@H]1O[C@H]([*])[C@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No beta-D-ribose sugar found"

    # Check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 2:
        return False, f"Insufficient hydroxyl groups ({hydroxyl_matches})"

    # Comprehensive nucleobase patterns
    base_patterns = [
        # Pyrimidines
        "[nX3]1c(=O)[nH]c(=O)cc1",  # Uracil
        "[nX3]1c(=O)[nH]c(=S)cc1",  # Thiouracil
        "[nX3]1c(N)nc(=O)cc1",      # Cytosine
        "[nX3]1c(=O)nc(N)cc1",      # Cytosine tautomer
        "[nX3]1cc[nH]c(=O)1",       # Basic pyrimidine
        
        # Purines
        "c1nc2c([nH]1)nc[nH]c2=O",  # Guanine
        "c1nc2c([nH]1)ncnc2N",      # Adenine
        "c1nc2c([nH]1)[nH]c(=O)nc2=O",  # Xanthine
        "c1[nH]c2ncnc-2n1",         # Basic purine
        
        # Modified bases
        "[nX3]1c(=O)[nH]c(=O)c(F)c1",  # Fluoropyrimidine
        "[nX3]1c(=O)[nH]c(=O)c(Br)c1", # Bromopyrimidine
        "[nX3]1c(=O)[nH]c(=O)c(C)c1",  # Methylpyrimidine
        "[nX3]1c(=O)[nH]c(=O)c(CC)c1", # Modified pyrimidine
        "c1nc2c([nH]1)nc(N)nc2=O",     # Modified purine
    ]

    has_base = False
    for pattern in base_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_base = True
            break
            
    if not has_base:
        return False, "No recognized nucleobase found"

    # Check for proper N-glycosidic bond with correct stereochemistry
    n_glycosidic_pattern = Chem.MolFromSmarts("[#7][C@H]1[O][C@H]([CH2][OH])[C@H]([OH])[C@@H]1[OH]")
    if not mol.HasSubstructMatch(n_glycosidic_pattern):
        return False, "No proper N-glycosidic bond with correct stereochemistry"

    # Ring count check
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient ring count"

    # Success case
    return True, "Contains beta-D-ribose sugar connected to nucleobase via N-glycosidic bond"