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

    # Exclude patterns (things we don't want to match)
    exclude_patterns = [
        # Carboxylic acids
        "[CX3](=O)[OX2H1]",
        # Esters
        "[CX3](=O)[OX2][CX4]",
        # Hemiacetals/hemiketal patterns (common in sugars)
        "[OX2H1][CX4]1[OX2][CX4][CX4][CX4][CX4]1",
        # Alpha-keto acids
        "[CX3](=O)[CX3](=O)[OX2H1]"
    ]
    
    for pattern in exclude_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            return False, "Contains excluded functional group pattern"

    # Core patterns for alpha-hydroxy ketone
    patterns = [
        # Basic pattern with explicit valence and connectivity
        "[OX2H1]-[CX4;H1;!$(C(O)(O))]-[CX3](=[OX1])-[#6;!$(C=O)]",
        
        # Pattern for cyclic ketones
        "[OX2H1]-[CX4;H1](-[#6,#7,#8,#16])-[CX3]1(=[OX1])-[#6](-[#6,#7,#8,#16])@[#6,#7,#8,#16]@[#6,#7,#8,#16]-1",
        
        # Pattern for conjugated ketones
        "[OX2H1]-[CX4;H1](-[#6,#7,#8,#16])-[CX3](=[OX1])-[CX3]=[CX3]",
        
        # Pattern for alpha-hydroxy ketones in complex ring systems
        "[OX2H1]-[CX4;H1;R](-[#6,#7,#8,#16])-[CX3;R](=[OX1])-[@#6]"
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
                
                # Validate alpha carbon
                if alpha_c.GetTotalNumHs() != 1:
                    continue
                    
                # Count non-H neighbors of alpha carbon
                non_h_neighbors = [n for n in alpha_c.GetNeighbors() 
                                 if n.GetAtomicNum() != 1]
                if len(non_h_neighbors) != 3:  # Must have exactly 3 non-H neighbors
                    continue

                # Verify one neighbor is OH, one is C=O, and one is C/heteroatom
                neighbor_types = set()
                for n in non_h_neighbors:
                    if n.GetAtomicNum() == 8 and n.GetTotalNumHs() == 1:  # OH
                        neighbor_types.add('OH')
                    elif n.GetAtomicNum() == 6 and any(nb.GetAtomicNum() == 8 and nb.GetIsAromatic() == False 
                                                      for nb in n.GetNeighbors()):  # C=O
                        neighbor_types.add('C=O')
                    else:  # Other group
                        neighbor_types.add('R')
                
                if neighbor_types != {'OH', 'C=O', 'R'}:
                    continue

                # Additional check for conjugated systems
                if any(a.GetIsAromatic() for a in [hydroxy_o, alpha_c, carbonyl_c]):
                    continue

                return True, "Contains secondary alpha-hydroxy ketone group"

    return False, "No secondary alpha-hydroxy ketone pattern found"