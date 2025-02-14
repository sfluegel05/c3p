"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: CHEBI:52526 alkanesulfonate oxoanion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule contains a sulfonate group (-SO3-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a sulfonate group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sulfonate group pattern
    sulfonate_pattern = Chem.MolFromSmarts("[S+3(=O)(-[O-])(-[O-])]")
    sulfonate_matches = mol.GetSubstructMatches(sulfonate_pattern)
    
    if len(sulfonate_matches) == 0:
        return False, "No sulfonate group found"
    
    # Check for at least one sulfonate group
    return True, f"Molecule contains {len(sulfonate_matches)} sulfonate group(s)"