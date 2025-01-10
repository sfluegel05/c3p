"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a fatty acid containing only methyl branches.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,OX2-]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carboxylic acid groups - should have just one
    carboxylic_matches = len(mol.GetSubstructMatches(carboxylic_pattern))
    if carboxylic_matches != 1:
        return False, f"Invalid number of carboxylic acid groups ({carboxylic_matches})"
    
    # Check for aromatic character
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"

    # Count carbons and check for valid atoms
    atom_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        
    if any(symbol not in ['C', 'H', 'O', 'N'] for symbol in atom_counts.keys()):
        return False, "Contains unexpected heteroatoms"

    # Identify branching patterns
    # First, find all carbons that have more than 2 carbon neighbors
    branched_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
            if len(carbon_neighbors) > 2:
                branched_carbons.append(atom.GetIdx())

    # Check each branch to ensure it's a methyl group
    for branch_point in branched_carbons:
        atom = mol.GetAtomWithIdx(branch_point)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon neighbor
                # Count carbons in this branch
                visited = set([atom.GetIdx()])
                branch_carbons = []
                stack = [(neighbor, [neighbor.GetIdx()])]
                
                while stack:
                    current, path = stack.pop()
                    if current.GetIdx() not in visited:
                        visited.add(current.GetIdx())
                        branch_carbons.append(current.GetIdx())
                        
                        # Don't follow path through the branch point
                        for next_atom in current.GetNeighbors():
                            if next_atom.GetIdx() not in visited and next_atom.GetIdx() != branch_point:
                                stack.append((next_atom, path + [next_atom.GetIdx()]))
                
                # If branch has more than 1 carbon and isn't part of main chain
                if len(branch_carbons) > 1 and not any(idx in branch_carbons for idx in mol.GetSubstructMatches(carboxylic_pattern)[0]):
                    return False, "Contains non-methyl branches"

    # Verify presence of at least one methyl branch
    methyl_pattern = Chem.MolFromSmarts("[CH3][CX4,CR4]")  # Methyl attached to sp3 or sp2 carbon
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "No methyl branches found"

    # Find main chain length (approximate using molecular formula)
    c_count = atom_counts.get('C', 0)
    if c_count < 4:
        return False, "Carbon chain too short"
    if c_count > 30:
        return False, "Carbon chain too long for typical fatty acid"

    # Check for disqualifying features
    disqualifying_patterns = [
        "[N+]",  # Quaternary nitrogen
        "[S,P,B,Si]",  # Other heteroatoms
        "C(=O)OC(=O)",  # Anhydrides
        "[OH]C(=O)[OH]"  # Carbonic acid
    ]
    
    for pattern in disqualifying_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, f"Contains disqualifying structural feature: {pattern}"

    return True, "Contains carboxylic acid group with methyl branching"