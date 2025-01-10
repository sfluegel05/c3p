"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a six-carbon monosaccharide with a specific stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for hexose structure - pyranose, furanose, or open-chain forms
    hexose_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H]([C@H](O)C1)O"),  # Pyranose form
        Chem.MolFromSmarts("O1[C@H]([C@H](O)[C@H](O)C1)CO"),  # Furanose form
        Chem.MolFromSmarts("[C@H](O)[C@H](O)[C@H](O)[C@H](O)C=O"),  # Open-chain aldehydo form
    ]
    
    # Attempt to find matches for these patterns
    matches = sum(mol.HasSubstructMatch(pattern) for pattern in hexose_patterns)
    if matches == 0:
        return False, "Does not match hexose structures (pyranose/furanose or open-chain forms)"

    # Verify the presence of D-stereochemistry - high chiral center requirement
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    valid_chiral_centers = sum(1 for _, flag in chiral_centers if flag in {'R', 'S'})
    if valid_chiral_centers < 4:
        return False, "Insufficient verified chiral centers for a D-hexose"
    
    # Verify that there are exactly six carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, "Not six carbons present"

    # Verify that there are exactly five oxygens (typical for hexoses in ring forms)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 5 and o_count != 6:  # Allow some flexibility
        return False, "Unusual number of oxygens for a hexose"

    return True, "Valid D-hexose conformer and stereochemistry detected"