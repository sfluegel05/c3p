"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: secondary alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has a carbonyl group adjacent to a carbon bearing
    one hydroxy group, one hydrogen, and one organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Multiple SMARTS patterns to catch different arrangements
    patterns = [
        # Basic alpha-hydroxy ketone pattern
        "[OX2H1,OX2-]-[CX4;H1](-[#6])-[CX3](=[OX1])",
        # Alternative pattern with different connectivity
        "[CX3](=[OX1])-[CX4;H1](-[OX2H1,OX2-])-[#6]",
        # Pattern for cyclic variants
        "[OX2H1,OX2-]-[CX4;H1](-[#6]@[#6])-[CX3](=[OX1])",
        # Pattern catching ring-based alpha-hydroxy ketones
        "[CX3]1(=[OX1])-[CX4;H1](-[OX2H1,OX2-])-[#6]-[#6]-[#6]-1"
    ]

    for pattern in patterns:
        pat = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(pat):
            matches = mol.GetSubstructMatches(pat)
            for match in matches:
                # Get the matched atoms
                hydroxy_o = mol.GetAtomWithIdx(match[0])
                alpha_c = mol.GetAtomWithIdx(match[1])
                carbonyl_c = mol.GetAtomWithIdx(match[3])

                # Verify basic requirements
                if alpha_c.GetTotalValence() != 4:  # Must be sp3
                    continue
                    
                # Count non-H neighbors of alpha carbon
                non_h_neighbors = sum(1 for neighbor in alpha_c.GetNeighbors() 
                                   if neighbor.GetAtomicNum() != 1)
                if non_h_neighbors != 3:  # Must have exactly 3 non-H neighbors
                    continue

                # Verify carbonyl is really a ketone (not aldehyde, acid, etc)
                carbonyl_neighbors = [atom for atom in carbonyl_c.GetNeighbors() 
                                   if atom.GetAtomicNum() != 8]
                if len(carbonyl_neighbors) != 1:
                    continue

                # Check if alpha carbon has one H (implicitly or explicitly)
                if alpha_c.GetTotalNumHs() != 1:
                    continue

                return True, "Contains secondary alpha-hydroxy ketone group"

    return False, "No secondary alpha-hydroxy ketone pattern found"