"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: CHEBI:35518 alkanesulfonate oxoanion

An alkanesulfonate oxoanion is defined as an alkanesulfonate in which the carbon at position 1
is attached to R, which can represent hydrogens, a carbon chain, or other groups.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for alkanesulfonate group pattern (S(=O)(=O)(=O)-C)
    alkanesulfonate_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])(=[O])-[C]")
    if not mol.HasSubstructMatch(alkanesulfonate_pattern):
        return False, "No alkanesulfonate group found"

    # Check if the molecule is an oxoanion
    formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if formal_charge != -1:
        return False, "Molecule is not an oxoanion (formal charge is not -1)"

    return True, "Molecule is an alkanesulfonate oxoanion"