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
    
    # Look for carbonyl carbon with one H and one alkyl/aryl group attached
    aldehyde_pattern = Chem.MolFromSmarts("[C;!$(C=O)]=[O;H]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    if not aldehyde_matches:
        return False, "No aldehyde functional group found"
    
    # Check for presence of at least one carbon
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 1:
        return False, "No carbon atoms found"
    
    # Check for exactly one aldehyde group
    if len(aldehyde_matches) > 1:
        return False, f"Found {len(aldehyde_matches)} aldehyde groups, need exactly 1"
    
    # Exclude carboxylic acids, esters, and amides
    excluded_patterns = [
        Chem.MolFromSmarts("C(=O)[O;H,-]"), # Carboxylic acid
        Chem.MolFromSmarts("C(=O)[O;!H0]"), # Ester
        Chem.MolFromSmarts("C(=O)[N]")      # Amide
    ]
    
    for pattern in excluded_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Molecule contains a non-aldehyde carbonyl group"
    
    # Check for presence of at least one alkyl/aryl group
    alkyl_pattern = Chem.MolFromSmarts("[C;!$(C=O)]")
    alkyl_matches = mol.GetSubstructMatches(alkyl_pattern)
    if not alkyl_matches:
        return False, "No alkyl/aryl group found"
    
    return True, "Contains an aldehyde functional group (R-C(=O)H)"