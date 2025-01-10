"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phenyl group (benzene ring)
    phenyl_pattern = Chem.MolFromSmarts('c1ccccc1')
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "No phenyl group (benzene ring) found"

    # Look for a three-carbon chain connected to a phenyl group
    propenyl_pattern = Chem.MolFromSmarts('c-C=C-C')
    if not mol.HasSubstructMatch(propenyl_pattern):
        return False, "No phenylpropanoid-like 3-carbon chain found"

    # Count number of oxygens and hydroxyl groups, typical in phenylpropanoids
    hydroxyl_pattern = Chem.MolFromSmarts('[OH]')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if len(hydroxyl_matches) == 0 and n_oxygens == 0:
        return False, "No functional groups typical of phenylpropanoids found"

    return True, "Has phenyl group, 3-carbon chain, and typical functional groups"