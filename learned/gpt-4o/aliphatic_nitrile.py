"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile contains a nitrile group derived from an aliphatic compound.

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

    # Look for nitrile group (C#N)
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    if not mol.HasSubstructMatch(nitrile_pattern):
        return False, "No nitrile group found"

    # Check if nitrile group's carbon is sp3 hybridized or part of an aliphatic chain
    for match in mol.GetSubstructMatches(nitrile_pattern):
        cn_atom, n_atom = match # cn_atom is the carbon in C#N
        carbon_atom = mol.GetAtomWithIdx(cn_atom)

        # Verify that the carbon is NOT directly bound to an sp2 or aromatic atom
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetHybridization() == Chem.HybridizationType.SP2 or neighbor.GetIsAromatic():
                return False, "Nitrile group is attached to a non-aliphatic group"

    return True, "Nitrile group is part of an aliphatic chain"