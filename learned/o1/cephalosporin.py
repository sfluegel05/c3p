"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: CHEBI:23066 cephalosporin
"""

from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins are beta-lactam antibiotics characterized by a beta-lactam ring fused to a six-membered dihydrothiazine ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define cephalosporin core SMARTS pattern
    cepha_core_smarts = '[#6]1[#6]2N1C(=O)[#6]=[#6][#16][#6]2'
    cepha_core_mol = Chem.MolFromSmarts(cepha_core_smarts)
    if cepha_core_mol is None:
        return False, "Failed to create cephalosporin core SMARTS pattern"

    # Check for cephalosporin core in the molecule
    if not mol.HasSubstructMatch(cepha_core_mol):
        return False, "Cephalosporin core structure not found"

    return True, "Contains cephalosporin core structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:23066',
        'name': 'cephalosporin',
        'definition': 'A class of beta-lactam antibiotics differing from the penicillins in having a 6-membered, rather than a 5-membered, side ring.',
    },
}