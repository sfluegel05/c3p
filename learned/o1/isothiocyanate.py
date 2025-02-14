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
    where R is any organic group (excluding hydrogen), and the nitrogen is bonded
    via a single bond to R, and double bonded to C in N=C=S.
    
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
    
    # Define isothiocyanate SMARTS pattern (nitrogen double-bonded to carbon, which is double-bonded to sulfur)
    isothiocyanate_pattern = Chem.MolFromSmarts("N=C=S")
    if isothiocyanate_pattern is None:
        return False, "Invalid SMARTS pattern"
    
    # Find all N=C=S groups
    matches = mol.GetSubstructMatches(isothiocyanate_pattern)
    if not matches:
        return False, "Does not contain isothiocyanate functional group (N=C=S)"
    
    # Check each match to ensure nitrogen is attached to an appropriate R group
    for match in matches:
        n_idx, c_idx, s_idx = match  # Indices of nitrogen, carbon, sulfur atoms
        n_atom = mol.GetAtomWithIdx(n_idx)
        c_atom = mol.GetAtomWithIdx(c_idx)
        s_atom = mol.GetAtomWithIdx(s_idx)
        
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
        
        # R group should be carbon
        if r_atom.GetAtomicNum() != 6:
            continue  # R group is not carbon, skip
        
        # Check that the bond between N and R group is a single bond
        bond = mol.GetBondBetweenAtoms(n_idx, r_atom.GetIdx())
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue  # Bond is not single, skip
        
        # Exclude cases where R group is a carbonyl carbon (C=O)
        is_carbonyl = False
        for bond in r_atom.GetBonds():
            neighbor = bond.GetOtherAtom(r_atom)
            if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                is_carbonyl = True  # R atom is a carbonyl carbon
                break
        if is_carbonyl:
            continue  # Skip if R atom is carbonyl carbon
        
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