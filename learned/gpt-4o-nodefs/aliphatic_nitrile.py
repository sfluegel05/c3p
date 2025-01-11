"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile contains a nitrile group (C#N) connected to an aliphatic carbon chain.

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

    # Look for nitrile group pattern (-C#N)
    nitrile_pattern = Chem.MolFromSmarts("[$([CX2]#N)]")
    matches = mol.GetSubstructMatches(nitrile_pattern)
    
    if not matches:
        return False, "No nitrile group found"

    # Check that nitrile group is not part of an aromatic system
    for match in matches:
        atom = mol.GetAtomWithIdx(match[0])
        # Check if the atom is part of an aliphatic chain
        if not atom.IsInRing() and atom.GetIsAromatic() == False:
            return True, "Contains an aliphatic nitrile group"

    return False, "Nitrile group is part of an aromatic system"