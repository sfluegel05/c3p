"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    Polyamines are characterized by the presence of multiple spaced amine groups across the molecule.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure presence of multiple amine groups
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;+0,!$([N]C=O)]")
    amine_matches = mol.GetSubstructMatches(amine_pattern)

    if len(amine_matches) < 2:
        return False, f"Found {len(amine_matches)} amine group(s), need at least 2"
    
    # Check paths for sufficient dispersion of amines over a carbon chain
    amine_atom_idxs = [idx for match in amine_matches for idx in match]
    significant_distances = []

    for i in range(len(amine_atom_idxs)):
        for j in range(i + 1, len(amine_atom_idxs)):
            path_length = Chem.rdmolops.GetShortestPath(mol, amine_atom_idxs[i], amine_atom_idxs[j])
            carbon_count = sum(1 for atom_idx in path_length if mol.GetAtomWithIdx(atom_idx).GetSymbol() == "C")
            
            # Considering significant spacing if a longer carbon chain (e.g., 3 carbons)
            if carbon_count >= 3:
                significant_distances.append((i, j, carbon_count))
    
    if not significant_distances:
        return False, "Amine groups are not sufficiently spaced by carbon chains"

    return True, "Contains multiple spaced amine groups with significant dispersion"