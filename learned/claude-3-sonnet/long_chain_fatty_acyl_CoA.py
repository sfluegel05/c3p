"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def get_longest_chain(mol, start_idx):
    """Find the longest continuous carbon chain from starting atom"""
    def dfs(atom_idx, visited=None):
        if visited is None:
            visited = set()
        
        if atom_idx in visited:
            return []
        
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        
        # Only continue if it's carbon and has proper connectivity
        if atom.GetAtomicNum() != 6 or atom.IsInRing():
            return []
        
        longest_path = [atom_idx]
        max_length = 0
        
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 16:  # Stop at sulfur
                continue
            if neighbor.GetAtomicNum() != 6:  # Skip non-carbon atoms
                continue
            
            path = dfs(neighbor.GetIdx(), visited.copy())
            if len(path) > max_length:
                max_length = len(path)
                longest_path = [atom_idx] + path
                
        return longest_path
    
    return dfs(start_idx)

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    Long-chain fatty acyl-CoAs have:
    - CoA moiety
    - Thioester linkage
    - Linear fatty acid chain length C13-C22
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define and validate all SMARTS patterns
    patterns = {
        'adenine': 'n1cnc2c(N)ncnc12',
        'thioester': 'C(=O)SC',
        'pantetheine': 'CCNC(=O)CCNC(=O)',
        'diphosphate': 'OP(=O)(O)OP(=O)(O)',
        'sugar': '[C,O]1[C,O][C,O][C,O][C,O]1',  # Generic sugar ring pattern
        'complex_ring': '[C,O]1[C,O][C,O][C,O][C,O][C,O]1',  # 6-membered rings
        'carboxyl': 'C(=O)[OH]'
    }
    
    smarts_patterns = {}
    for name, pattern in patterns.items():
        smarts_mol = Chem.MolFromSmarts(pattern)
        if smarts_mol is None:
            return False, f"Internal error: Invalid SMARTS pattern for {name}"
        smarts_patterns[name] = smarts_mol

    # Check for required CoA structural features
    required_features = ['adenine', 'thioester', 'pantetheine', 'diphosphate']
    for feature in required_features:
        if not mol.HasSubstructMatch(smarts_patterns[feature]):
            return False, f"Missing {feature} moiety required for CoA structure"

    # Check for disqualifying features
    if mol.HasSubstructMatch(smarts_patterns['sugar']):
        return False, "Contains sugar moiety - not a simple fatty acyl-CoA"
    
    if mol.HasSubstructMatch(smarts_patterns['complex_ring']):
        return False, "Contains complex ring structure - not a simple fatty acyl-CoA"
    
    if len(mol.GetSubstructMatches(smarts_patterns['carboxyl'])) > 0:
        return False, "Contains carboxyl group - not a simple fatty acyl-CoA"

    # Find the thioester carbon
    thioester_matches = mol.GetSubstructMatches(smarts_patterns['thioester'])
    if not thioester_matches:
        return False, "Could not identify thioester linkage"
    
    # Get the carbon atom connected to the thioester sulfur
    thioester_carbon = thioester_matches[0][0]
    
    # Get the main carbon chain
    chain_atoms = get_longest_chain(mol, thioester_carbon)
    chain_length = len(chain_atoms)
    
    # Check chain length (C13-C22)
    if chain_length < 13:
        return False, f"Fatty acid chain too short (C{chain_length}, need C13-C22)"
    if chain_length > 22:
        return False, f"Fatty acid chain too long (C{chain_length}, need C13-C22)"

    # Count double bonds in the chain
    double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in chain_atoms and end_idx in chain_atoms:
                double_bonds += 1

    # Look for modifications
    hydroxy_pattern = Chem.MolFromSmarts('CO')
    oxo_pattern = Chem.MolFromSmarts('CC(=O)C')
    
    features = []
    if double_bonds > 0:
        features.append(f"{double_bonds} double bond(s)")
    if mol.HasSubstructMatch(hydroxy_pattern):
        features.append("hydroxy group(s)")
    if mol.HasSubstructMatch(oxo_pattern):
        features.append("oxo group(s)")
    
    feature_str = ", ".join(features)
    if feature_str:
        feature_str = f" with {feature_str}"
    
    return True, f"Long-chain fatty acyl-CoA (C{chain_length}){feature_str}"