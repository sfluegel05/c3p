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

    # Get the main carbon skeleton
    main_skeleton = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts('[!#6]'))
    c_count = main_skeleton.GetNumAtoms()
    
    # Hexose must have exactly 6 carbons in the main skeleton
    if c_count != 6:
        return False, f"Found {c_count} carbons in main skeleton, need exactly 6"

    # Check for aldehyde group (C=O at position 1 in linear form)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CX4]")
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains aldehyde group at position 1 (aldohexose)"

    # Check for ketone group (C=O at position 2 in linear form)
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # For cyclic forms, check if the molecule has characteristic hexose patterns
    # More general pattern for cyclic hexoses
    cyclic_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@H](O)[C@H](O)[C@H](O1)CO)")
    if mol.HasSubstructMatch(cyclic_pattern):
        # Check if it has a potential ketone group
        for match in ketone_matches:
            # Check if the ketone is on a carbon that could be position 2
            atom = mol.GetAtomWithIdx(match[0])
            neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
            if 8 in neighbors:  # Check if connected to oxygen (carbonyl)
                return True, "Cyclic hexose form with potential ketone group detected"
        return True, "Cyclic hexose form detected"

    # Check for ketone in linear form
    if ketone_matches:
        # Try to find a ketone at position 2 in a linear chain
        for match in ketone_matches:
            atom = mol.GetAtomWithIdx(match[0])
            # Check if it's in a linear chain and could be position 2
            if len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]) == 2:
                return True, "Contains ketone group at position 2 (ketohexose)"

    # If neither aldehyde nor ketone is found, it's not a hexose
    return False, "No aldehyde or ketone group found in linear or cyclic form"