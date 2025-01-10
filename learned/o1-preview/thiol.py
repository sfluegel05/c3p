"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound in which a thiol group (-SH) is attached 
    to a carbon atom of any aliphatic or aromatic moiety.

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

    # Add explicit hydrogens to ensure accurate matching
    mol = Chem.AddHs(mol)

    # Define SMARTS pattern for thiol group (-SH attached to carbon)
    thiol_smarts = "[#6]-[SH]"  # Carbon atom bonded to sulfur with a hydrogen
    thiol_pattern = Chem.MolFromSmarts(thiol_smarts)

    # Search for thiol group in the molecule
    if mol.HasSubstructMatch(thiol_pattern):
        return True, "Contains thiol group (-SH) attached to carbon atom"

    return False, "No thiol group (-SH) attached to carbon atom found"