"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: CHEBI:51766 long-chain fatty alcohol
A fatty alcohol with a chain length ranging from C13 to C22.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Lipinski import LengthDescriptor

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a chain length ranging from C13 to C22 with an -OH group attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for an -OH group
    oh_pattern = Chem.MolFromSmarts("[OX1H]")
    oh_indices = mol.GetSubstructMatches(oh_pattern)
    if not oh_indices:
        return False, "No hydroxy (-OH) group found"

    # Count length of longest aliphatic chain
    length_descriptor = LengthDescriptor(LengthDescriptor.LMCHANGESLD)
    longest_chain_length = length_descriptor(mol)
    if not (13 <= longest_chain_length <= 22):
        return False, f"Longest aliphatic chain length is {longest_chain_length}, not in the range C13 to C22"

    # Check if -OH group is attached to the longest aliphatic chain
    for oh_idx in oh_indices:
        oh_atom = mol.GetAtomWithIdx(oh_idx)
        if any(bond.GetOtherAtomIdx(oh_idx) in length_descriptor.GetAtomIndices(mol) for bond in oh_atom.GetBonds()):
            break
    else:
        return False, "Hydroxy (-OH) group not attached to the longest aliphatic chain"

    # Check for disallowed functional groups
    disallowed_patterns = [
        Chem.MolFromSmarts("[C$(C=O)O]"),  # Esters
        Chem.MolFromSmarts("C(=O)C"),  # Ketones
        Chem.MolFromSmarts("C(=O)O"),  # Carboxylic acids
    ]
    for pattern in disallowed_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains disallowed functional groups"

    return True, "Molecule contains a straight-chain aliphatic alcohol with a chain length between C13 and C22"