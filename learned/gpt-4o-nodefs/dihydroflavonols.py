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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for dihydroflavanone core with chiral centers
    dihydroflavonol_core_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](Oc2cc(O)cc(O)c2C1=O)c1ccccc1")  # Possible dihydroflavanone core
    if not mol.HasSubstructMatch(dihydroflavonol_core_pattern):
        return False, "No dihydroflavonol core found"
    
    # Check for necessary chiral centers
    stereo = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(stereo) < 2:
        return False, "Chiral centers do not match typical dihydroflavonol structure"
    
    # Look for multiple hydroxyl groups - basic pattern coverage
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    n_hydroxyls = len(hydroxyl_matches)
    if n_hydroxyls < 3:
        return False, f"Found {n_hydroxyls} hydroxyl groups, expected at least 3"
    
    return True, "Contains dihydroflavonol features including core structure, chiral centers, and appropriate hydroxylation"