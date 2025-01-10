"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion contains a sulfonate group (SO3-) attached to an alkane,
    which implies non-participation in olefinic or aromatic bonds directly at the attachment point.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES into a molecular object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improve SMARTS patterns to specify context around the sulfonate
    # The carbon directly bonded to S should not be involved in double bonds (ensure no ketones, etc.)
    carbon_sulfonate_pattern = Chem.MolFromSmarts("[CX4][S]([O-])(=O)=O")
    if mol.HasSubstructMatch(carbon_sulfonate_pattern):
        return True, "Contains a sulfonate group appropriately bonded to an alkane-like carbon"

    return False, "No matching alkanesulfonate oxoanion structure found"

# Example SMILES for testing
smiles_list = [
    "[O-]S(C[C@H]([C@@H](O)C=O)O)(=O)=O",
    "CCS([O-])(=O)=O",
    "CCCCS([O-])(=O)=O",
    # Adding some corner cases to ensure robustness
    "C(S([O-])(=O)=O)=C",       # Should return False due to olefinic carbon
    "c1ccccc1S([O-])(=O)=O",    # Should return False due to aromatic involvement
]

for smiles in smiles_list:
    result, reason = is_alkanesulfonate_oxoanion(smiles)
    print(f"SMILES: {smiles} -> Result: {result}, Reason: {reason}")