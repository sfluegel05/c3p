"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile contains a nitrile group derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for nitrile group (C#N)
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)

    if len(nitrile_matches) == 0:
        return False, "No nitrile group found"

    for match in nitrile_matches:
        cn_atom, n_atom = match  # cn_atom is the carbon in C#N
        carbon_atom = mol.GetAtomWithIdx(cn_atom)
        
        # Initialize a flag for aliphatic nature
        is_aliphatic = True

        # Check neighbors of the nitrile carbon for aromaticity/sp2 hybridization
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIdx() == n_atom:  # Skip the nitrogen atom in the nitrile
                continue
            if neighbor.GetIsAromatic() or neighbor.GetHybridization() != Chem.HybridizationType.SP3:
                is_aliphatic = False
                break

        # Check further bonds for potentially non-aliphatic characteristics
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIdx() == n_atom or not is_aliphatic:
                continue
            if any(next_neighbor.GetIsAromatic() or next_neighbor.GetHybridization() != Chem.HybridizationType.SP3
                    for next_neighbor in neighbor.GetNeighbors() if next_neighbor.GetIdx() != cn_atom):
                is_aliphatic = False
                break

        # Return result based on aliphatic check
        if is_aliphatic:
            return True, "Nitrile group is part of an aliphatic chain"

    return False, "Nitrile group is attached to a non-aliphatic group"