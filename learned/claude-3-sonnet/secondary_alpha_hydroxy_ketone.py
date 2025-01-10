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

    # Core patterns for alpha-hydroxy ketone
    patterns = [
        # Basic pattern: OH-CH(R)-C(=O)-R
        "[OX2H1]-[CX4;H1](-[#6,#7,#8,#16])-[CX3](=[OX1])-[#6]",
        
        # Reversed pattern: R-C(=O)-CH(OH)-R
        "[#6]-[CX3](=[OX1])-[CX4;H1](-[OX2H1])-[#6,#7,#8,#16]",
        
        # Cyclic pattern where ketone is in ring
        "[OX2H1]-[CX4;H1](-[#6,#7,#8,#16])-[CX3]1(=[OX1])-[#6,#7,#8,#16]@[#6,#7,#8,#16]@[#6,#7,#8,#16]-1",
        
        # Pattern where OH is deprotonated
        "[OX2-]-[CX4;H1](-[#6,#7,#8,#16])-[CX3](=[OX1])-[#6]",
        
        # Pattern for ring junction cases
        "[OX2H1]-[CX4;H1](-[#6,#7,#8,#16]@[#6,#7,#8,#16])-[CX3](=[OX1])-[@#6,#7,#8,#16]"
    ]

    for pattern in patterns:
        pat = Chem.MolFromSmarts(pattern)
        if not pat:
            continue
            
        if mol.HasSubstructMatch(pat):
            matches = mol.GetSubstructMatches(pat)
            for match in matches:
                # Get the matched atoms
                hydroxy_o = mol.GetAtomWithIdx(match[0])
                alpha_c = mol.GetAtomWithIdx(match[1])
                carbonyl_c = mol.GetAtomWithIdx(match[2])
                
                # Basic validation of alpha carbon
                if alpha_c.GetTotalNumHs() != 1:
                    continue
                    
                # Count non-H neighbors of alpha carbon
                non_h_neighbors = sum(1 for neighbor in alpha_c.GetNeighbors() 
                                   if neighbor.GetAtomicNum() != 1)
                if non_h_neighbors < 3:  # Must have at least 3 non-H neighbors
                    continue

                # Verify carbonyl is ketone-like (not carboxylic acid, ester, etc)
                carbonyl_neighbors = carbonyl_c.GetNeighbors()
                has_valid_ketone = False
                for neighbor in carbonyl_neighbors:
                    if neighbor.GetAtomicNum() == 6:  # Carbon
                        non_o_neighbors = sum(1 for n in neighbor.GetNeighbors() 
                                           if n.GetAtomicNum() != 8)
                        if non_o_neighbors >= 1:
                            has_valid_ketone = True
                            break
                
                if not has_valid_ketone:
                    continue

                return True, "Contains secondary alpha-hydroxy ketone group"

    return False, "No secondary alpha-hydroxy ketone pattern found"