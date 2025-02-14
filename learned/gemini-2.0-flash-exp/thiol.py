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

    # Look for thiol group (-SH)
    thiol_pattern = Chem.MolFromSmarts("[SH1]")
    if not mol.HasSubstructMatch(thiol_pattern):
       return False, "No thiol (-SH) group found"

    # Check if -SH is attached to a carbon
    carbon_thiol_pattern = Chem.MolFromSmarts("[CX4,cx3][SH1]")
    if not mol.HasSubstructMatch(carbon_thiol_pattern):
        return False, "Thiol group not directly attached to a carbon"

    return True, "Molecule contains a thiol group attached to a carbon"