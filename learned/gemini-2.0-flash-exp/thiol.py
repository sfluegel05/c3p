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

    # Look for a sulfur atom attached to a carbon atom (aliphatic or aromatic)
    thiol_pattern_aliphatic = Chem.MolFromSmarts("[S][CX4]")
    thiol_pattern_aromatic = Chem.MolFromSmarts("[S][c]")

    aliphatic_matches = mol.GetSubstructMatches(thiol_pattern_aliphatic)
    aromatic_matches = mol.GetSubstructMatches(thiol_pattern_aromatic)

    if not aliphatic_matches and not aromatic_matches:
          return False, "No sulfur atom attached to a carbon atom found"


    for match in aliphatic_matches:
        sulfur_atom = mol.GetAtomWithIdx(match[0])
        if sulfur_atom.GetTotalNumHs() <= 1:
             return True, "Molecule contains a thiol group attached to a carbon"

    for match in aromatic_matches:
        sulfur_atom = mol.GetAtomWithIdx(match[0])
        if sulfur_atom.GetTotalNumHs() <= 1:
             return True, "Molecule contains a thiol group attached to a carbon"
             
    return False, "No thiol group (-SH) attached to a carbon atom found"