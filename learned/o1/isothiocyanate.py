"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: isothiocyanate
"""
from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    An isothiocyanate is an organosulfur compound with the general formula R-N=C=S,
    where R is any organic group (excluding hydrogen).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define isothiocyanate SMARTS pattern
    # Nitrogen double-bonded to carbon, which is double-bonded to sulfur
    isothiocyanate_pattern = Chem.MolFromSmarts("N=C=S")
    if isothiocyanate_pattern is None:
        return False, "Invalid SMARTS pattern"
    
    # Find all N=C=S groups
    matches = mol.GetSubstructMatches(isothiocyanate_pattern)
    if not matches:
        return False, "Does not contain isothiocyanate functional group (N=C=S)"
    
    # Check each match to ensure nitrogen is attached to an R group (non-hydrogen)
    for match in matches:
        n_idx, c_idx, s_idx = match  # Indices of nitrogen, carbon, sulfur atoms
        n_atom = mol.GetAtomWithIdx(n_idx)
        c_atom = mol.GetAtomWithIdx(c_idx)
        s_atom = mol.GetAtomWithIdx(s_idx)
        
        # Check that nitrogen is not in a ring
        if n_atom.IsInRing():
            continue  # Skip if nitrogen is in a ring

        # Nitrogen should have degree 2 (connected to carbon in N=C=S and R group)
        if n_atom.GetDegree() != 2:
            continue  # Skip if nitrogen does not have exactly two neighbors

        # Get the neighbor atom of nitrogen that is not the carbon in N=C=S
        neighbor_atoms = [atom for atom in n_atom.GetNeighbors() if atom.GetIdx() != c_idx]
        if len(neighbor_atoms) != 1:
            continue  # Nitrogen does not have exactly one R group attached

        r_atom = neighbor_atoms[0]
        # Check that R group is not hydrogen
        if r_atom.GetAtomicNum() == 1:
            continue  # R group is hydrogen, skip

        # Check that the bond between N and R group is a single bond
        bond = mol.GetBondBetweenAtoms(n_idx, r_atom.GetIdx())
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue  # Bond is not single, skip

        # Passed all checks, is an isothiocyanate
        return True, "Contains isothiocyanate functional group (R-N=C=S)"
    
    return False, "Does not contain isothiocyanate functional group (R-N=C=S)"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:47905',
        'name': 'isothiocyanate',
        'definition': 'An organosulfur compound with the general formula R-N=C=S.'
    }
}