"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: CHEBI:16456 Aldehyde
An aldehyde is a compound RC(=O)H, in which a carbonyl group is bonded to one hydrogen atom and to one R group.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carbonyl carbon with one H and one R group attached
    aldehyde_pattern = Chem.MolFromSmarts("C(=O)[H]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    if not aldehyde_matches:
        return False, "No aldehyde functional group found"
    
    # Check for presence of at least one carbon
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 1:
        return False, "No carbon atoms found"
    
    # Check for presence of only one aldehyde group
    if len(aldehyde_matches) > 1:
        return False, f"Found {len(aldehyde_matches)} aldehyde groups, need exactly 1"
    
    # Check if the carbonyl carbon is not part of a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Molecule contains a carboxylic acid group, not an aldehyde"
    
    return True, "Contains an aldehyde functional group (R-C(=O)H)"