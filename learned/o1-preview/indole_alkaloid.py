"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: indole alkaloid
"""

from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid is an alkaloid containing an indole skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define indole SMARTS pattern (matches indole and N-substituted indoles)
    indole_smarts = 'c1c[nH,n]c2cccc2c1'
    indole = Chem.MolFromSmarts(indole_smarts)

    # Check if molecule has indole substructure
    if not mol.HasSubstructMatch(indole):
        return False, "No indole skeleton found"

    # Get nitrogen atoms in indole substructure(s)
    indole_matches = mol.GetSubstructMatches(indole)
    indole_nitrogen_indices = set()
    for match in indole_matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7:
                indole_nitrogen_indices.add(idx)

    # Get all nitrogen atoms in the molecule
    nitrogen_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    total_nitrogens = len(nitrogen_indices)
    indole_nitrogens = len(indole_nitrogen_indices)
    additional_nitrogens = total_nitrogens - indole_nitrogens

    if additional_nitrogens >= 1:
        return True, "Contains indole skeleton and additional nitrogen atoms"
    else:
        return False, "Contains indole skeleton but no additional nitrogen atoms"

__metadata__ = {'chemical_class': { 'id': '',  # CHEBI ID can be added if known
                             'name': 'indole alkaloid',
                             'definition': 'An alkaloid containing an indole skeleton.',
                             'parents': []},
    # Additional metadata can be added here if needed
}