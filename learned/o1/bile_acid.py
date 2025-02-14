"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: Bile Acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid is a hydroxy-5beta-cholanic acid with a steroid backbone,
    hydroxyl groups at specific positions, and a carboxylic acid side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens (explicit)
    mol = Chem.AddHs(mol)
    
    # Define the 5beta-cholanic acid skeleton using SMILES
    # The SMILES represents the steroid nucleus with the 5beta-configuration
    cholanic_acid_smiles = "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@H]3CC[C@@H]4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    skeleton_pattern = Chem.MolFromSmiles(cholanic_acid_smiles)
    if skeleton_pattern is None:
        return False, "Error in skeleton pattern"

    # Match the steroid skeleton with correct stereochemistry
    if not mol.HasSubstructMatch(skeleton_pattern, useChirality=True):
        return False, "Molecule does not match 5beta-cholanic acid skeleton"

    # Check for at least one hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    num_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if num_hydroxyls == 0:
        return False, "No hydroxyl groups found"

    # Check for carboxylic acid group (COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    return True, "Molecule matches bile acid structural features"