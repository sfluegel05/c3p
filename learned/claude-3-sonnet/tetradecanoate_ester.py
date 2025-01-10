"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: tetradecanoate ester
A fatty acid ester obtained by condensation of the carboxy group of tetradecanoic acid 
with a hydroxy group of an alcohol or phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule contains a tetradecanoate ester group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for ester group pattern
    ester_pattern = Chem.MolFromSmarts("[CX4,CH2X4:1][CX3:2](=[OX1:3])[OX2:4][C,c:5]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if not ester_matches:
        return False, "No ester group found"
    
    # For each ester group, check if it's part of a tetradecanoyl group
    for match in ester_matches:
        # Get the carbon atom of the ester carbonyl
        carbonyl_carbon = mol.GetAtomWithIdx(match[1])
        
        # Count carbons in the chain before the ester group
        visited = set()
        chain = []
        current = mol.GetAtomWithIdx(match[0])  # Start from first carbon
        
        while current is not None and len(chain) < 14:
            if current.GetIdx() in visited:
                break
                
            visited.add(current.GetIdx())
            if current.GetAtomicNum() == 6:  # Carbon
                chain.append(current)
                
            # Get next carbon in chain
            next_atom = None
            for neighbor in current.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    next_atom = neighbor
                    break
            current = next_atom
            
        # Check if we found exactly 13 carbons in the chain (plus the carbonyl carbon = 14)
        if len(chain) == 13:
            # Verify it's a straight chain
            is_straight = True
            for atom in chain:
                if len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]) > 2:
                    is_straight = False
                    break
                    
            if is_straight:
                return True, "Contains tetradecanoate (myristoyl) ester group"
    
    return False, "No tetradecanoate ester group found"