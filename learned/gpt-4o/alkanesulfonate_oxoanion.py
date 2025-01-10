"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion contains a sulfonate group (SO3-) attached to an alkane
    or a carbon-containing group.

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

    # Define improved SMARTS patterns
    # This pattern detects a carbon attached directly to a sulfonate oxoanion structure
    carbon_sulfonate_pattern = Chem.MolFromSmarts("[CX4,CX3,CX2,CX1,CX0][S]([O-])(=O)=O")

    # Check for the pattern matches
    if mol.HasSubstructMatch(carbon_sulfonate_pattern):
        return True, "Contains a sulfonate group appropriately bonded to a carbon"

    return False, "No matching sulfonate structure found indicative of an oxoanion"

# Example SMILES for testing
smiles_list = [
    "[O-]S(C[C@H]([C@@H](O)C=O)O)(=O)=O",
    "CCS([O-])(=O)=O",
    "CCCCS([O-])(=O)=O",
]

for smiles in smiles_list:
    result, reason = is_alkanesulfonate_oxoanion(smiles)
    print(f"SMILES: {smiles} -> Result: {result}, Reason: {reason}")