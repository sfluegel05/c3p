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

    # Look for structure typical of phenylpropanoids, like 3-carbon chain
    # connected to phenyl group with possible modifications
    three_carbon_chain_connected = Chem.MolFromSmarts('c-[C!H0]-[C!H0]-[C!H0]')
    if not mol.HasSubstructMatch(three_carbon_chain_connected):
        return False, "No 3-carbon chain connected to phenyl group found"

    # Count the number of oxygens or typical functional groups (e.g., hydroxyl, carbonyl)
    relevant_oxygen_substructures = [
        Chem.MolFromSmarts('[OH]'),   # Hydroxyl
        Chem.MolFromSmarts('[CX3]=[OX1]'),  # Carbonyl in ketones/aldehydes
        Chem.MolFromSmarts('[CX3](=O)[OX2H1]'),  # Carboxyl group
        Chem.MolFromSmarts('[OX2][CX4]')  # Ether linkages
    ]

    oxygen_functions = sum(
        mol.HasSubstructMatch(functional_group) for functional_group in relevant_oxygen_substructures
    )
    
    if oxygen_functions == 0:
        return False, "No typical functional groups found"

    return True, "Contains phenyl group, 3-carbon chain, and functional groups typical of phenylpropanoids"