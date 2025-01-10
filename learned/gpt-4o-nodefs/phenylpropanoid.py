"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem

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
        return None, "Invalid SMILES string"

    # Check for phenyl group (benzene ring)
    phenyl_group = Chem.MolFromSmarts('c1ccccc1')
    if not mol.HasSubstructMatch(phenyl_group):
        return False, "No phenyl group (benzene ring) found"

    # Look for a 3-carbon chain attached to phenyl group
    propanoid_chain = Chem.MolFromSmarts('c1c[c,C][C!H0][C!H0]c1')
    if not mol.HasSubstructMatch(propanoid_chain):
        return False, "No appropriate 3-carbon chain connected to phenyl group found"

    # Check for common functional groups such as hydroxyl, carbonyl
    functional_groups = [
        Chem.MolFromSmarts('[OH]'),   # Hydroxyl
        Chem.MolFromSmarts('[CX3]=[OX1]'),  # Carbonyl in ketones/aldehydes
        Chem.MolFromSmarts('[CX3](=O)[OX2H1]'),  # Carboxyl group
        Chem.MolFromSmarts('[OX2][CX4]'),  # Ether linkages
        Chem.MolFromSmarts('[CX4][OH]'),   # Alcohols
        Chem.MolFromSmarts('[OX2][CX3]=[CX3]') # Propenoic acid esters
    ]

    if not any(mol.HasSubstructMatch(fg) for fg in functional_groups):
        return False, "No characteristic functional groups found"

    return True, "Contains phenyl group, 3-carbon chain, and functional groups typical of phenylpropanoids"