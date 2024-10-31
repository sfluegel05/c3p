from rdkit import Chem
from rdkit.Chem import AllChem

def is_uric_acid(smiles: str):
    """
    Determines if a molecule is a uric acid derivative (oxopurine that is an oxidation product of purine metabolism)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a uric acid derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for basic purine core (fused 5-6 rings with nitrogens)
    purine_pattern = Chem.MolFromSmarts("[#7]1[#6]2[#7][#6][#7][#6][#6]2[#7][#6]1")
    if not mol.HasSubstructMatch(purine_pattern):
        purine_pattern2 = Chem.MolFromSmarts("[#7]1[#6]2[#7][#6][#7][#6][#7]2[#6][#6]1")
        if not mol.HasSubstructMatch(purine_pattern2):
            return False, "Does not contain purine core structure"

    # Count oxo groups (=O) and hydroxy groups (-OH) in various positions
    patterns = [
        "[O,OH][#6]1[#7][#6]2[#7][#6][#7][#6][#6,#7]2[#7][#6]1", # Pattern for 2-position
        "[#7]1[#6]2[#7][#6]([O,OH])[#7][#6][#6,#7]2[#7][#6]1",   # Pattern for 6-position
        "[#7]1[#6]2[#7][#6][#7][#6]([O,OH])[#6,#7]2[#7][#6]1"    # Pattern for 8-position
    ]

    oxidized_positions = set()
    for i, pattern in enumerate(patterns):
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        if matches:
            oxidized_positions.add(i)

    if len(oxidized_positions) < 3:
        return False, f"Not fully oxidized - found oxidation at positions {list(oxidized_positions)}, need all three positions"

    # Count nitrogen atoms to ensure correct structure
    n_pattern = Chem.MolFromSmarts("[#7]")
    n_count = len(mol.GetSubstructMatches(n_pattern))
    if n_count < 4:
        return False, "Insufficient nitrogen atoms for purine structure"

    return True, f"Valid uric acid with oxidation at positions 2,6,8"
# Pr=None
# Recall=0.0