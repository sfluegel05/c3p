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

    # First check for ester group
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][C,c]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Find all ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    for match in ester_matches:
        carbonyl_carbon = mol.GetAtomWithIdx(match[0])
        
        # Now check if this carbonyl carbon has a 13-carbon chain attached
        # (making it a C14 chain in total)
        chain_length = 0
        visited = set()
        current_atom = carbonyl_carbon
        
        # Traverse the carbon chain
        while True:
            visited.add(current_atom.GetIdx())
            
            # Find next carbon in chain
            carbon_neighbors = [n for n in current_atom.GetNeighbors() 
                             if n.GetAtomicNum() == 6 and 
                             n.GetIdx() not in visited and
                             not n.IsInRing()]  # Exclude ring carbons
            
            if len(carbon_neighbors) != 1:
                break
                
            current_atom = carbon_neighbors[0]
            chain_length += 1
            
            # Check for branching
            if len([n for n in current_atom.GetNeighbors() 
                   if n.GetAtomicNum() == 6]) > 2:
                break
        
        if chain_length == 13:  # We found a 13-carbon chain attached to carbonyl
            # Verify it's saturated (no double/triple bonds)
            is_saturated = True
            for atom_idx in visited:
                atom = mol.GetAtomWithIdx(atom_idx)
                if any(bond.GetBondType() != Chem.BondType.SINGLE 
                      for bond in atom.GetBonds()):
                    is_saturated = False
                    break
            
            if is_saturated:
                return True, "Contains tetradecanoate (myristoyl) ester group"
    
    return False, "No tetradecanoate ester group found"