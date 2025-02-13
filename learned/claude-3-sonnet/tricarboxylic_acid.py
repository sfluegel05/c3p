"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: CHEBI:33569 tricarboxylic acid
An oxoacid containing three carboxy groups.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid should have three carboxyl groups (-COOH) attached
    to a common carbon backbone or ring system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxyl group pattern (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("[C](=O)(O)[OH]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) != 3:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups, need exactly 3"

    # Find the carbon atoms attached to the carboxyl groups
    carboxyl_carbons = set(mol.GetAtomWithIdx(match[0]).GetNeighbors()[0].GetIdx() for match in carboxyl_matches)

    # Check if the carboxyl groups are connected to a common carbon backbone or ring
    backbone_atoms = set()
    for carbon in carboxyl_carbons:
        backbone_atoms.update(Chem.AtomAtomPathFinder.findAtomPathWalk(mol, carbon, maxPath=3))

    # If all carboxyl groups are connected, the intersection of their paths should be non-empty
    if len(backbone_atoms.intersection(carboxyl_carbons)) >= 3:
        return True, "Contains three carboxyl groups attached to a common carbon backbone or ring"
    else:
        return False, "Carboxyl groups are not connected to a common carbon backbone or ring"