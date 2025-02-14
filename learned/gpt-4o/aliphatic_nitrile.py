"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is any nitrile derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for nitrile group (-C#N)
    nitrile_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
    if not mol.HasSubstructMatch(nitrile_pattern):
        return False, "No nitrile group found"

    for match in mol.GetSubstructMatches(nitrile_pattern):
        carbon_idx, nitrogen_idx = match
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)

        # Determine aliphatic property by examining all neighbors of the carbon in the nitrile group
        is_aliphatic = True
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                is_aliphatic = False
                break

        if is_aliphatic:
            return True, "Nitrile group is derived from an aliphatic compound"

    # If no aliphatic nitrile groups are found
    return False, "Nitrile group is not derived from an aliphatic compound"