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

    # Check that nitrile group is part of an aliphatic context
    for match in mol.GetSubstructMatches(nitrile_pattern):
        carbon_idx, nitrogen_idx = match
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)

        # Ensure carbon atom is not directly aromatic
        if carbon_atom.GetIsAromatic():
            return False, "Nitrile group carbon is aromatic"

        # Analyze the chain linked to nitrile carbon
        aliphatic_chain = False
        for neighbor in carbon_atom.GetNeighbors():
            # Check if neighbor is part of an aliphatic chain
            if not neighbor.GetIsAromatic():
                aliphatic_chain = True
                break

        if aliphatic_chain:
            return True, "Nitrile group is derived from an aliphatic compound"

    # If none of the nitrile linkages are truly aliphatic
    return False, "Nitrile group is not aliphatic"