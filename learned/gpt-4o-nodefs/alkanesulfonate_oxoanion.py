"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion contains a sulfonate group attached to an aliphatic carbon chain.

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

    # Define sulfonate group pattern
    sulfonate_pattern = Chem.MolFromSmarts("S([O-])(=O)=O")
    if not mol.HasSubstructMatch(sulfonate_pattern):
        return False, "Missing sulfonate group S([O-])(=O)=O"

    # Check for aliphatic carbon connection
    aliphatic_carbon_pattern = Chem.MolFromSmarts("[CX4]")
    for match in mol.GetSubstructMatches(sulfonate_pattern):
        sulfur_atom_idx = match[0]  # Index of Sulfur atom in the match
        sulfur_atom = mol.GetAtomWithIdx(sulfur_atom_idx)

        # Iterate over neighbors of the sulfur atom
        for neighbor in sulfur_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and mol.GetAtomWithIdx(neighbor.GetIdx()).GetIsAromatic() is False:
                # Check if the neighbor carbon atom is an aliphatic one
                if mol.HasSubstructMatch(aliphatic_carbon_pattern, atoms=[neighbor.GetIdx()]):
                    return True, "Sulfonate group attached to an aliphatic structure"
    
    return False, "Sulfonate group not connected to an aliphatic structure"