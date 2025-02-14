"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: hydroxy fatty acid
"""

from rdkit import Chem

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is any fatty acid carrying one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid group found"

    # Get indices of hydroxyl oxygen atoms in carboxylic acid
    carboxy_oxygen_indices = [match[2] for match in carboxy_matches]

    # Find all hydroxy groups (excluding carboxylic acid hydroxyl)
    hydroxy_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)

    # Remove hydroxy groups that are part of carboxylic acid
    hydroxy_indices = [match[0] for match in hydroxy_matches if match[0] not in carboxy_oxygen_indices]

    if not hydroxy_indices:
        return False, "No hydroxy substituents found in the molecule excluding carboxylic acid group"

    # Check for aliphatic chain connected to carboxylic acid carbon
    carboxylic_carbon_indices = [match[0] for match in carboxy_matches]
    has_aliphatic_chain = False
    for c_idx in carboxylic_carbon_indices:
        carbon_atom = mol.GetAtomWithIdx(c_idx)
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != c_idx:
                has_aliphatic_chain = True
                break
        if has_aliphatic_chain:
            break

    if not has_aliphatic_chain:
        return False, "No aliphatic chain connected to carboxylic acid group"

    return True, "Molecule is a hydroxy fatty acid"