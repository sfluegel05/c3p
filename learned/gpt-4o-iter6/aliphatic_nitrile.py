"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile contains a nitrile (-C#N) group attached to an aliphatic chain.

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

    # Look for nitrile group pattern (C#N)
    nitrile_pattern = Chem.MolFromSmarts("[#6]#N")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    if not nitrile_matches:
        return False, "Nitrile group (-C#N) not found"
    
    # Check if nitrile group is attached to an aliphatic carbon
    for match in nitrile_matches:
        # Get the carbon atom in the C#N group
        carbon_index = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_index)
        # Check if the carbon is attached directly to an sp3 carbon (aliphatic)
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                return True, "Nitrile group attached to aliphatic carbon"
    
    return False, "Nitrile group not attached to aliphatic chain"