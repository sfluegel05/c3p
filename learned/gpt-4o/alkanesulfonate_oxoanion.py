"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion contains a sulfonate group (SO3-) attached to an alkane or a carbon-containing group.

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
    
    # Define the sulfonate group pattern connected to carbon
    sulfonate_pattern = Chem.MolFromSmarts("CS([O-])(=O)=O")
    
    # Check if the molecule has the sulfonate group
    if mol.HasSubstructMatch(sulfonate_pattern):
        return True, "Contains a sulfonate group bonded to a carbon atom"
    else:
        # Check the pattern "[R]S([O-])(=O)=O" where R is a carbon atom
        carbon_sulfonate_pattern = Chem.MolFromSmarts("[CX4,CX3]S([O-])(=O)=O")
        if mol.HasSubstructMatch(carbon_sulfonate_pattern):
            return True, "Contains a sulfonate group bonded to a carbon atom"
        else:
            return False, "No sulfonate group attached to a carbon found"

# Example SMILES for testing
smiles_list = [
    "[O-]S(C[C@H]([C@@H](O)C=O)O)(=O)=O",
    "CCS([O-])(=O)=O",
    "CCCCS([O-])(=O)=O",
]

for smiles in smiles_list:
    result, reason = is_alkanesulfonate_oxoanion(smiles)
    print(f"SMILES: {smiles} -> Result: {result}, Reason: {reason}")