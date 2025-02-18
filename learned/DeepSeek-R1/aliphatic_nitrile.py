"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: CHEBI:18311 aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile has nitrile groups (-C#N) attached to aliphatic carbons (non-aromatic).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    # Parse SMILES and sanitize to ensure proper aromaticity detection
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all nitrile groups (C#N) with adjacent carbon context
    # Pattern matches: [C(connected to other atoms)]#[N(terminal)]-[any atom]
    nitrile_pattern = Chem.MolFromSmarts('[CX3]#[NX1]-[*]')
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    
    if not nitrile_matches:
        return False, "No nitrile groups found"
    
    # Check the carbon adjacent to each nitrile group
    for match in nitrile_matches:
        # Match format: (nitrile_carbon_idx, nitrile_nitrogen_idx, adjacent_atom_idx)
        adjacent_atom = mol.GetAtomWithIdx(match[2])
        
        # Check if adjacent atom is aromatic (part of aromatic system)
        if adjacent_atom.GetIsAromatic():
            return False, f"Nitrile group attached to aromatic system at position {adjacent_atom.GetIdx()+1}"

    return True, "All nitrile groups are attached to aliphatic carbons"