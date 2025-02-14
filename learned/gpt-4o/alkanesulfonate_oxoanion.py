"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for alkanesulfonate (-SO3^-) structure: S([O-])(=O)=O
    sulfonate_pattern = Chem.MolFromSmarts("S([O-])(=O)=O")
    if not mol.HasSubstructMatch(sulfonate_pattern):
        return False, "No sulfonate group S([O-])(=O)=O found"

    # Check if the sulfonate group is attached to a carbon (i.e., an alkane chain)
    # This checks if the sulfur atom (attached to S([O-])(=O)=O) has a further carbon 'R' group.
    for match in mol.GetSubstructMatches(sulfonate_pattern):
        sulfur_atom = mol.GetAtomWithIdx(match[0])
        # Check connected atoms to sulfur for at least one connected carbon
        connected_carbons = [n.GetAtomicNum() == 6 for n in sulfur_atom.GetNeighbors()]
        if any(connected_carbons):
            return True, "Contains a sulfonate group attached to a carbon chain or other R group"

    return False, "Sulfonate group not attached to a carbon chain or proper R group"