"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: olefinic fatty acid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    An olefinic fatty acid is any fatty acid containing at least one C=C double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
        return False, f"Only {c_count} carbon atoms found, too few for a fatty acid"

    # Check for at least one C=C double bond
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) == 0:
        return False, "No C=C double bonds found"

    return True, "Contains carboxylic acid group and at least one C=C double bond"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'olefinic fatty acid',
        'definition': 'Any fatty acid containing at least one C=C double bond.',
        'parents': []
    }
}