"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile contains a nitrile (-C#N) group bonded to a carbon
    that is part of an aliphatic (non-aromatic, may be cyclic) structure.

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

    # Define nitrile pattern
    nitrile_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)

    if not nitrile_matches:
        return False, "Nitrile group (-C#N) not found"
    
    # Check if nitrile group is bonded to an aliphatic carbon
    for match in nitrile_matches:
        carbon_idx = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # Ensure the carbon atom is not aromatic
        if not carbon_atom.GetIsAromatic():
            # Look at the carbon atom's neighbors
            neighbors = carbon_atom.GetNeighbors()
            for neighbor in neighbors:
                # Check that the neighbor is a carbon and not aromatic
                if neighbor.GetAtomicNum() == 6 and not neighbor.GetIsAromatic():
                    # Check for appropriate hybridization in neighbors to confirm aliphatic structure
                    if (neighbor.GetHybridization() == Chem.rdchem.HybridizationType.SP3 or 
                        neighbor.GetHybridization() == Chem.rdchem.HybridizationType.SP2):
                        return True, "Nitrile group attached to non-aromatic aliphatic carbon"

    return False, "Nitrile group not appropriately part of an aliphatic chain"