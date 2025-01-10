"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    The bond between carbons 11 and 12 (counting from CoA end) must be saturated.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an 11,12-saturated fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES with stereochemistry
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Essential CoA structural elements
    # Adenine base
    adenine_pattern = Chem.MolFromSmarts('n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine base found"

    # Phosphate groups - look for multiple patterns
    phosphate_patterns = [
        Chem.MolFromSmarts('OP([O-])([O-])=O'),
        Chem.MolFromSmarts('OP([O-])(=O)O'),
        Chem.MolFromSmarts('P([O-])([O-])(=O)O')
    ]
    
    phosphate_count = sum(len(mol.GetSubstructMatches(pattern)) 
                         for pattern in phosphate_patterns)
    if phosphate_count < 3:
        return False, f"Insufficient phosphate groups for CoA, found {phosphate_count}"

    # Check for negative charges
    negative_charges = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
    if negative_charges < 4:
        return False, f"Insufficient negative charges for CoA(4-), found {negative_charges}"

    # Look for thioester linkage - more flexible pattern
    thioester_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[SX2][CH2]')
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Look for pantetheine fragments - split into parts for more flexible matching
    amide_pattern = Chem.MolFromSmarts('NC(=O)CC')
    if len(mol.GetSubstructMatches(amide_pattern)) < 2:
        return False, "Missing pantetheine amide groups"

    # Find the thioester carbon and trace the fatty acid chain
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Could not locate thioester group"

    # Get the carbon attached to sulfur
    thioester_idx = thioester_matches[0][0]  # First carbon of thioester
    
    # Find the fatty acid chain by following carbons
    visited = set()
    chain = []
    current = thioester_idx
    
    while True:
        visited.add(current)
        chain.append(current)
        
        # Get carbon neighbors not yet visited
        atom = mol.GetAtomWithIdx(current)
        c_neighbors = [n.GetIdx() for n in atom.GetNeighbors() 
                      if n.GetAtomicNum() == 6 and n.GetIdx() not in visited]
        
        if not c_neighbors:
            break
            
        # Sort neighbors by chain length to follow main chain
        def get_chain_length(atom_idx):
            visited_local = visited.copy()
            count = 0
            current_local = atom_idx
            while True:
                visited_local.add(current_local)
                atom_local = mol.GetAtomWithIdx(current_local)
                c_neighbors_local = [n.GetIdx() for n in atom_local.GetNeighbors() 
                                   if n.GetAtomicNum() == 6 and n.GetIdx() not in visited_local]
                if not c_neighbors_local:
                    break
                current_local = c_neighbors_local[0]
                count += 1
            return count
        
        # Follow the longest chain
        current = max(c_neighbors, key=get_chain_length)
        
        # Check bond between positions 11 and 12
        if len(chain) == 11 and c_neighbors:
            bond_11_12 = mol.GetBondBetweenAtoms(chain[10], current)
            if bond_11_12.GetBondType() != Chem.BondType.SINGLE:
                return False, "Bond between carbons 11 and 12 is unsaturated"
    
    # Verify chain length
    if len(chain) < 12:
        return False, "Fatty acid chain too short (less than 12 carbons)"

    return True, "11,12 bond is saturated in fatty acyl-CoA chain"