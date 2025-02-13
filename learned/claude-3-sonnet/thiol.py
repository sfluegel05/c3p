"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound with a thiol group (-SH) attached to an
    aliphatic or aromatic carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a thiol, False otherwise
        str: Reason for classification

    Raises:
        ValueError: If the input SMILES string is invalid
    """

    try:
        mol = Chem.MolFromSmiles(smiles)
    except ValueError:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    # Find the sulfhydryl (-SH) group
    thiol_pattern = Chem.MolFromSmarts("[SH]")
    thiol_matches = mol.GetSubstructMatches(thiol_pattern)

    if not thiol_matches:
        return False, "No sulfhydryl (-SH) group found"

    # Check if the sulfur atom is connected to an aliphatic or aromatic carbon
    aliphatic_pattern = Chem.MolFromSmarts("[CH3,CH2,CH1]")
    aromatic_pattern = Chem.MolFromSmarts("[c]")

    for thiol_match in thiol_matches:
        sulfur_atom = mol.GetAtomWithIdx(thiol_match)
        neighbors = sulfur_atom.GetNeighbors()

        for neighbor in neighbors:
            if mol.GetAtomWithIdx(neighbor.GetIdx()).HasQuery(aliphatic_pattern) or \
               mol.GetAtomWithIdx(neighbor.GetIdx()).HasQuery(aromatic_pattern):
                return True, "Contains a sulfhydryl (-SH) group attached to an aliphatic or aromatic carbon"

    return False, "Sulfhydryl (-SH) group not attached to an aliphatic or aromatic carbon"