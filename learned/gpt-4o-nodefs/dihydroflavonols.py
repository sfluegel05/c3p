"""
Classifies: CHEBI:48039 dihydroflavonols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol typically contains a characteristic dihydroflavonone core structure
    with hydroxylation and specific chiral centers.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to obtain molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # More flexible SMARTS pattern for dihydroflavonol core allowing variations
    dihydroflavonol_core_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@@H](Oc2cc(O)cc(O)c2C1=O)c1ccccc1"),  # Basic core
        Chem.MolFromSmarts("O[C@H]1[C@H](Oc2cc(O)c(O)cc2C1=O)c1ccc(O)c(O)c1"), # Variant
    ]
    
    # Check for core structure matches
    core_match = any(mol.HasSubstructMatch(pattern) for pattern in dihydroflavonol_core_patterns)
    if not core_match:
        return False, "No dihydroflavonol core structure found"
    
    # Identify chiral centers and ensure they conform to typical dihydroflavonol chiral patterns
    stereo_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(stereo_centers) < 2:
        return False, "Insufficient chiral centers for dihydroflavonol config"
    
    # Check for necessary number of hydroxyl groups (usually more than three)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    n_hydroxyls = len(set(match[0] for match in hydroxyl_matches))
    if n_hydroxyls < 3:
        return False, f"Insufficient hydroxyl groups: expected at least 3, found {n_hydroxyls}"
    
    return True, "Contains dihydroflavonol features including core structure, chiral centers, and appropriate hydroxylation"