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
    nitrile_pattern = Chem.MolFromSmarts("[CX2]#N")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    
    if not nitrile_matches:
        return False, "No nitrile group found"

    # Check that nitrile group is connected to an aliphatic carbon (sp3, not part of any aromatic system)
    for match in nitrile_matches:
        nitrile_carbon = mol.GetAtomWithIdx(match[0])
        if nitrile_carbon.GetIsAromatic():
            continue

        # Consider the carbon directly bonded to the nitrile carbon
        for neighbor in nitrile_carbon.GetNeighbors():
            if neighbor.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and not neighbor.GetIsAromatic():
                return True, "Contains an aliphatic nitrile group"

    return False, "Nitrile group is not connected to an aliphatic chain"