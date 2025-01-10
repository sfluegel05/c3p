"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def count_chain_carbons(mol, start_atom_idx):
    """Helper function to count carbons in a chain from starting atom"""
    visited = set()
    def dfs(atom_idx, in_chain=True):
        if atom_idx in visited:
            return 0
        visited.add(atom_idx)
        
        atom = mol.GetAtomWithIdx(atom_idx)
        count = 1 if (in_chain and atom.GetAtomicNum() == 6) else 0
        
        for neighbor in atom.GetNeighbors():
            # Don't count carbons in CoA part
            if neighbor.GetAtomicNum() == 16:  # Sulfur
                continue
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited:
                count += dfs(neighbor_idx, in_chain)
        return count
    
    return dfs(start_atom_idx)

def count_double_bonds(mol, in_chain_atoms):
    """Helper function to count double bonds in the fatty acid chain"""
    count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Only count if both atoms are in the chain
            if bond.GetBeginAtomIdx() in in_chain_atoms and bond.GetEndAtomIdx() in in_chain_atoms:
                count += 1
    return count

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    Long-chain fatty acyl-CoAs have:
    - CoA moiety
    - Thioester linkage
    - Fatty acid chain length C13-C22
    
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
        'ribose_phosphate': 'OCC1OC(n)C(O)C1OP(O)(O)=O'
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

    # Find the thioester carbon
    thioester_matches = mol.GetSubstructMatches(smarts_patterns['thioester'])
    if not thioester_matches:
        return False, "Could not identify thioester linkage"
    
    # Get the carbon atom connected to the thioester sulfur
    thioester_carbon = thioester_matches[0][0]
    
    # Count carbons in the fatty acid chain
    chain_carbons = count_chain_carbons(mol, thioester_carbon)
    
    # Check chain length (C13-C22)
    if chain_carbons < 13:
        return False, f"Fatty acid chain too short (C{chain_carbons}, need C13-C22)"
    if chain_carbons > 22:
        return False, f"Fatty acid chain too long (C{chain_carbons}, need C13-C22)"

    # Get atoms in the chain for checking modifications
    visited = set()
    def get_chain_atoms(atom_idx):
        if atom_idx in visited:
            return set()
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        chain_atoms = {atom_idx}
        
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 16:  # Skip sulfur and beyond
                continue
            chain_atoms.update(get_chain_atoms(neighbor.GetIdx()))
        return chain_atoms
    
    chain_atoms = get_chain_atoms(thioester_carbon)
    
    # Count modifications
    double_bonds = count_double_bonds(mol, chain_atoms)
    
    # Look for common modifications
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
    
    return True, f"Long-chain fatty acyl-CoA (C{chain_carbons}){feature_str}"