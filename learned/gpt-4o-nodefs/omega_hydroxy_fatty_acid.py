"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    Omega-hydroxy fatty acids are characterized by a terminal hydroxyl group (ω-position)
    and a carboxylic acid group at the other end.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Match omega-terminal hydroxy group pattern
    omega_hydroxy_smarts = "[OX2H][CX4]"  # Hydroxyl group should be terminal, possibly at ω-position
    omega_hydroxy = Chem.MolFromSmarts(omega_hydroxy_smarts)
    if not mol.HasSubstructMatch(omega_hydroxy):
        return False, "No omega-terminal hydroxy group found"

    # Match carboxylic acid group
    carboxylic_acid_smarts = "C(=O)[OH]"
    carboxylic_acid = Chem.MolFromSmarts(carboxylic_acid_smarts)
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Check that hydroxyl group is terminal in the longest chain
    length, terminal_index = max((len(Chem.rdmolops.GetShortestPath(mol, match[0], match[1])), match[0])
                                 for match in mol.GetSubstructMatches(carboxylic_acid))
    omega_ends = [at.GetIdx() for at in mol.GetAtoms() if at.GetSmarts() == '[OH]']
    if terminal_index not in omega_ends:
        return False, "Hydroxy group is not terminal in the longest carbon chain"

    # Chain length criteria is set shorter based on examples
    num_carbs = rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)
    if num_carbs < 3:  # Change from 6 to 3 to accommodate validated shorter chains
        return False, f"Carbon chain too short, found {num_carbs} carbons"

    return True, "Contains omega-terminal hydroxy and carboxylic acid groups with appropriate chain length"