"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide with either an aldehyde group at position 1 (aldohexose)
    or a ketone group at position 2 (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons in the main skeleton
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Hexose must have exactly 6 carbons in the main skeleton
    if c_count != 6:
        return False, f"Found {c_count} carbons in main skeleton, need exactly 6"

    # Check for aldehyde group (C=O at position 1 in linear form)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CX4]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if aldehyde_matches:
        return True, "Contains aldehyde group at position 1 (aldohexose)"

    # Check for ketone group (C=O at position 2 in linear form)
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # For cyclic forms, check if the molecule can be linearized to show the ketone at position 2
    for match in ketone_matches:
        # Create a copy of the molecule
        mol_copy = Chem.Mol(mol)
        # Try to break rings to linearize the molecule
        Chem.Kekulize(mol_copy)
        # Check if the ketone is now at position 2
        if len(match) > 1 and match[0] != match[1]:
            return True, "Contains ketone group at position 2 (ketohexose)"

    # Check for cyclic forms that can open to reveal aldehyde/ketone
    # Look for characteristic hexose patterns
    hexose_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@H](O)[C@H](O)[C@H](O1)CO)")
    if mol.HasSubstructMatch(hexose_pattern):
        return True, "Cyclic hexose form detected"

    # If neither aldehyde nor ketone is found, it's not a hexose
    return False, "No aldehyde or ketone group found in linear or cyclic form"