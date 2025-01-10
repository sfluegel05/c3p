"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile contains a nitrile (-C#N) group bonded to a carbon that is part of an aliphatic (non-aromatic) structure.

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

    # Look for nitrile group pattern (C#N) bonded to non-aromatic carbons
    nitrile_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)

    if not nitrile_matches:
        return False, "Nitrile group (-C#N) not found"
    
    # Check if nitrile group is part of an aliphatic chain
    for match in nitrile_matches:
        carbon_index = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_index)
        # Ensure carbon is not aromatic and part of an aliphatic chain
        if not carbon_atom.GetIsAromatic():
            # Check neighbors to ensure part of an aliphatic chain or presence of sp3/sp2 hybridization
            for neighbor in carbon_atom.GetNeighbors():
                if neighbor.GetHybridization() in (Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP2):
                    return True, "Nitrile group attached to aliphatic carbon or part of an aliphatic chain"
    
    return False, "Nitrile group not part of an aliphatic chain"