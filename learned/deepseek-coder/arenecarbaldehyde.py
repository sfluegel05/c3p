"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: CHEBI:22680 arenecarbaldehyde
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is any aldehyde in which the carbonyl group is attached to an aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aldehyde group (-CHO)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group found"

    # Check if the aldehyde carbon is directly attached to an aromatic atom
    aromatic_aldehyde_pattern = Chem.MolFromSmarts("[a]-[CX3H1]=O")
    aromatic_aldehyde_matches = mol.GetSubstructMatches(aromatic_aldehyde_pattern)
    if not aromatic_aldehyde_matches:
        return False, "Aldehyde group not directly attached to an aromatic atom"

    # Check if the molecule contains at least one aromatic ring
    aromatic_ring_pattern = Chem.MolFromSmarts("[a]")
    aromatic_ring_matches = mol.GetSubstructMatches(aromatic_ring_pattern)
    if not aromatic_ring_matches:
        return False, "No aromatic ring found in the molecule"

    # Relax the restrictions on additional functional groups
    # Only exclude molecules where the additional functional groups interfere with the aldehyde-aromatic bond
    # For example, exclude molecules where the aldehyde is part of a larger functional group (e.g., carboxylic acid)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        # Check if the carboxylic acid is part of the aldehyde group
        for match in mol.GetSubstructMatches(carboxylic_acid_pattern):
            if match[0] in [m[0] for m in aldehyde_matches]:
                return False, "Aldehyde group is part of a carboxylic acid"

    return True, "Contains an aldehyde group directly attached to an aromatic atom"