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

    # Check for CoA pattern
    coa_pattern = Chem.MolFromSmarts('NC(=O)CCNC(=O)[CH]OC(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC([CH](O)[CH]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA group found"

    # Check for thioester linkage (-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts('C(=O)S')
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Find the thioester carbon and trace the fatty acid chain
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Could not locate thioester group"
        
    # Get the carbon atom index of the thioester C(=O)
    thioester_carbon = thioester_matches[0][0]
    
    # Build path from thioester carbon
    fatty_chain = []
    visited = set()
    current = thioester_carbon
    
    # Trace carbons in the fatty acid chain
    while True:
        visited.add(current)
        neighbors = [n for n in mol.GetAtomWithIdx(current).GetNeighbors() 
                    if n.GetAtomicNum() == 6 and n.GetIdx() not in visited]
        if not neighbors:
            break
        # Follow the longest carbon chain
        current = neighbors[0].GetIdx()
        fatty_chain.append(current)
        
        # Check bond between positions 11 and 12 if we have enough carbons
        if len(fatty_chain) >= 12:
            bond_11_12 = mol.GetBondBetweenAtoms(fatty_chain[10], fatty_chain[11])
            if bond_11_12.GetBondType() != Chem.BondType.SINGLE:
                return False, "Bond between carbons 11 and 12 is unsaturated"
    
    # Make sure chain is long enough
    if len(fatty_chain) < 12:
        return False, "Fatty acid chain too short (less than 12 carbons)"
        
    return True, "11,12 bond is saturated in fatty acyl-CoA chain"