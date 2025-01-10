"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:27208 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    VLCFAs have chain length >C22. Those >C27 are ultra-long-chain fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VLCFA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if carbon_count <= 22:
        return False, f"Chain length C{carbon_count} is too short (must be >C22)"
    
    # Find the longest carbon chain
    # First, get all carbons connected to the carboxyl group
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not matches:  # This should never happen as we checked above
        return False, "No carboxylic acid group found"
    
    carboxyl_carbon = matches[0][0]  # Get the carbon of the carboxyl group
    
    # Function to find longest chain from a starting atom
    def find_longest_chain(atom, visited=None):
        if visited is None:
            visited = set()
        
        visited.add(atom.GetIdx())
        max_length = 1
        
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:
                length = find_longest_chain(neighbor, visited.copy())
                max_length = max(max_length, length + 1)
                
        return max_length
    
    # Get the longest carbon chain length
    longest_chain = find_longest_chain(mol.GetAtomWithIdx(carboxyl_carbon))
    
    if longest_chain <= 22:
        return False, f"Longest continuous carbon chain (C{longest_chain}) is too short (must be >C22)"
    
    # Check for common VLCFA features
    is_ultra = longest_chain > 27
    has_double_bonds = mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]=[#6]"))
    has_triple_bonds = mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]#[#6]"))
    has_hydroxy = mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2H1][#6]"))
    has_cyclopropyl = mol.HasSubstructMatch(Chem.MolFromSmarts("[C]1[C][C]1"))
    has_branching = mol.HasSubstructMatch(Chem.MolFromSmarts("[#6][CH]([#6])[#6]"))
    
    # Build description
    features = []
    if is_ultra:
        features.append("ultra-long-chain")
    if has_double_bonds:
        features.append("unsaturated")
    if has_triple_bonds:
        features.append("contains triple bonds")
    if has_hydroxy:
        features.append("hydroxylated")
    if has_cyclopropyl:
        features.append("contains cyclopropyl rings")
    if has_branching:
        features.append("branched")
    
    feature_str = " and ".join(features) if features else "saturated"
    chain_desc = f"C{longest_chain}"
    
    return True, f"{chain_desc} {feature_str} fatty acid"