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

        # Check nitrile carbon atom environment
        # Disallow direct connection to aromatic systems
        if carbon_atom.GetIsAromatic():
            continue  # skip if carbon in the nitrile group is marked aromatic
        
        # Ensure surrounding atoms are not overly conjugated
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                continue  # skip this nitrile if immediately connected to aromatic
        
        # If it's non-aromatic or aliphatic, accept it
        return True, "Nitrile group is derived from an aliphatic compound"

    # If none of the nitrile contexts suggest aliphatic origins
    return False, "Nitrile group is not derived from an aliphatic compound"