"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: CHEBI:17478 aldehyde

Definition: A compound RC(=O)H, in which a carbonyl group is bonded to one hydrogen atom and to one R group.
"""

from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is characterized by a carbonyl group (C=O) bonded to one hydrogen atom and one R group.

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

    # Aldehyde pattern: carbonyl group with carbon bonded to hydrogen
    aldehyde_pattern = Chem.MolFromSmarts('[CH]=O')
    matches = mol.GetSubstructMatches(aldehyde_pattern)

    if matches:
        num_aldehyde_groups = len(matches)
        return True, f"Contains {num_aldehyde_groups} aldehyde functional group(s)"
    else:
        return False, "Does not contain aldehyde functional group"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17478',
        'name': 'aldehyde',
        'definition': 'A compound RC(=O)H, in which a carbonyl group is bonded to one hydrogen atom and to one R group.',
        'parents': []
    }
}