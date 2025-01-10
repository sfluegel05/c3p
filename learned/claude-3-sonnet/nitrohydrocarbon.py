"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: CHEBI:51753 nitrohydrocarbon
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon where one or more hydrogens are replaced by nitro groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of nitro groups (-NO2)
    nitro_pattern = Chem.MolFromSmarts('[$([N+](=[OX1])([O-])[#6])]')
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro groups found"

    # Get all atoms
    atoms = mol.GetAtoms()
    
    # Check if molecule contains only C, H, N, and O
    allowed_atoms = {1, 6, 7, 8}  # H, C, N, O atomic numbers
    for atom in atoms:
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Contains elements other than C, H, N, O"

    # Check for carbon atoms (must have at least one)
    carbon_atoms = [atom for atom in atoms if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return False, "No carbon atoms found"

    # Check for connected carbon framework
    carbon_smarts = Chem.MolFromSmarts('C~C')
    if mol.HasSubstructMatch(carbon_smarts):
        connected_carbons = True
    else:
        # Single carbon with nitro group is also valid
        if len(carbon_atoms) == 1 and len(nitro_matches) >= 1:
            connected_carbons = True
        else:
            connected_carbons = False
    
    if not connected_carbons:
        return False, "No connected carbon framework found"

    # Check for unwanted functional groups that would make it not a pure nitrohydrocarbon
    unwanted_groups = [
        ('C(=O)O', 'carboxylic acid'),
        ('CO', 'alcohol'),
        ('C(=O)', 'carbonyl'),
        ('CN', 'amine'),
        ('CS', 'thiol/sulfide'),
        ('CCl', 'chloride'),
        ('CBr', 'bromide'),
        ('CF', 'fluoride'),
        ('CI', 'iodide'),
        ('C[OH]', 'hydroxyl'),
    ]
    
    for pattern, group_name in unwanted_groups:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, f"Contains {group_name} group"

    # Count nitro groups
    num_nitro = len(nitro_matches)
    
    return True, f"Hydrocarbon with {num_nitro} nitro group{'s' if num_nitro > 1 else ''}"