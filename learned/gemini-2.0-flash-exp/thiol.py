"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound in which a thiol group (-SH) is attached to a carbon atom of any aliphatic or aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a sulfur atom with one single bond.
    thiol_pattern = Chem.MolFromSmarts("[S]")
    sulfur_atoms = mol.GetSubstructMatches(thiol_pattern)
    if not sulfur_atoms:
        return False, "No sulfur atom found"

    # Check if the sulfur atom is bonded to a carbon atom
    carbon_thiol_pattern = Chem.MolFromSmarts("[S][C]")
    if not mol.HasSubstructMatch(carbon_thiol_pattern):
        return False, "Thiol group not attached to a carbon atom"

    return True, "Molecule contains a thiol group attached to a carbon"