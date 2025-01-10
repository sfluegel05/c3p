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
        cn_atom, n_atom = match # cn_atom is the carbon in C#N
        carbon_atom = mol.GetAtomWithIdx(cn_atom)

        # Check if the nitrile carbon is only connected to other non-aromatic, non-sp2 carbon atoms
        is_aliphatic = True
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIdx() == n_atom:  # Skip the nitrogen atom in the nitrile
                continue
            if neighbor.GetIsAromatic() or neighbor.GetHybridization() == Chem.HybridizationType.SP2:
                is_aliphatic = False
                break
            if neighbor.GetHybridization() == Chem.HybridizationType.SP:
                # Check if the neighboring SP hybridized atom has only aliphatic connections
                for sub_neighbor in neighbor.GetNeighbors():
                    if sub_neighbor.GetIdx() != cn_atom and (
                        sub_neighbor.GetIsAromatic() or sub_neighbor.GetHybridization() == Chem.HybridizationType.SP2):
                        is_aliphatic = False
                        break
        
        if is_aliphatic:
            return True, "Nitrile group is part of an aliphatic chain"

    return False, "Nitrile group is attached to a non-aliphatic group"