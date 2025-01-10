"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    A monoterpenoid indole alkaloid contains an indole moiety derived from L-tryptophan and
    a monoterpene-derived moiety (C10 unit), usually linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for indole moiety (allowing substitutions at nitrogen)
    indole_pattern = Chem.MolFromSmarts('c1cccc2c1nccc2')  # Indole core with substituted N
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole moiety found"

    # Check for monoterpene-derived moiety (approximate)
    # Monoterpenes are C10 units built from isoprene units (C5H8)
    # We can check if the molecule has at least 20 carbons (C10 from indole + C10 from monoterpene)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 20:
        return False, f"Too few carbons ({num_carbons}) to be a monoterpenoid indole alkaloid"
        
    # Check for connectivity between indole and monoterpene units
    # For simplicity, we'll accept presence of indole and sufficient size

    return True, "Contains indole moiety and sufficient size to be a monoterpenoid indole alkaloid"

__metadata__ = {
    'chemical_class': {
        'name': 'monoterpenoid indole alkaloid',
        'definition': 'A terpenoid indole alkaloid which is biosynthesised from L-tryptophan and diisoprenoid (usually secologanin) building blocks.'
    },
    'config': {
        # Configuration parameters if any
    }
}