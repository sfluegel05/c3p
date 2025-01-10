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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for key CoA structural elements
    # Adenine base with ribose
    adenine_pattern = Chem.MolFromSmarts('n1cnc2c(N)ncnc12')
    # Phosphates with negative charges
    phosphate_pattern = Chem.MolFromSmarts('[O-]P([O-])(=O)O')
    # Pantetheine part with thioester
    pantetheine_pattern = Chem.MolFromSmarts('NC(=O)CCNC(=O)CCSC(=O)')

    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine base found"
    if len(mol.GetSubstructMatches(phosphate_pattern)) < 3:
        return False, "Missing required phosphate groups"
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No pantetheine-thioester linkage found"

    # Find the thioester carbon and trace the fatty acid chain
    thioester_pattern = Chem.MolFromSmarts('SC(=O)C')
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Could not locate thioester group"

    # Get the carbon atoms in the fatty acid chain
    start_carbon = thioester_matches[0][2]  # Carbon attached to C=O
    
    # Trace the main carbon chain
    chain = []
    current = start_carbon
    visited = set()
    
    while True:
        visited.add(current)
        chain.append(current)
        
        # Get carbon neighbors not yet visited
        c_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(current).GetNeighbors() 
                      if n.GetAtomicNum() == 6 and n.GetIdx() not in visited]
        
        if not c_neighbors:
            break
            
        # Follow the first carbon neighbor
        current = c_neighbors[0]
        
        # Check bond between positions 11 and 12 when we reach them
        if len(chain) == 11 and c_neighbors:
            bond_11_12 = mol.GetBondBetweenAtoms(chain[10], c_neighbors[0])
            if bond_11_12.GetBondType() != Chem.BondType.SINGLE:
                return False, "Bond between carbons 11 and 12 is unsaturated"
    
    # Verify chain length
    if len(chain) < 12:
        return False, "Fatty acid chain too short (less than 12 carbons)"

    return True, "11,12 bond is saturated in fatty acyl-CoA chain"