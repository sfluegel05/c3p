"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile contains a nitrile (-C#N) group bonded to a carbon
    that is part of an aliphatic (non-aromatic, non-aromatic rings included) structure.

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

    # Define nitrile pattern: carbon triple-bonded to nitrogen
    nitrile_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)

    if not nitrile_matches:
        return False, "Nitrile group (-C#N) not found"
    
    # Check if nitrile group is bonded to an aliphatic carbon
    for match in nitrile_matches:
        # The carbon in the nitrile triple bond is indexed at 0
        c_nitrile_idx = match[0]
        c_nitrile_atom = mol.GetAtomWithIdx(c_nitrile_idx)
        
        # Check neighbors of the nitrile carbon
        neighbors = c_nitrile_atom.GetNeighbors()
        for neighbor in neighbors:
            # Ensure the neighbor is a carbon
            if neighbor.GetAtomicNum() == 6:
                # Ensure this neighbor carbon is not aromatic
                if not neighbor.GetIsAromatic():
                    # Check if the bonding structure resembles an aliphatic chain
                    if all(not atom.GetIsAromatic() for atom in neighbor.GetNeighbors()):
                        return True, "Nitrile group attached to an aliphatic carbon"
                    else:
                        return False, "Neighbor carbons part of aromatic structure"
    
    return False, "Nitrile group not adequately part of an aliphatic chain"