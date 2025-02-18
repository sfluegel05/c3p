"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
"""
Classifies: Straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid.
    Criteria:
    - Contains a carboxylic acid group
    - All carbon-carbon bonds are single (saturated)
    - No branching in the carbon chain
    - No additional carbon substituents (only straight chain)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for carboxylic acid group (-COOH)
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group"
    
    # Check all carbon-carbon bonds are single (saturated)
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        if begin.GetSymbol() == 'C' and end.GetSymbol() == 'C' and bond.GetBondType() != Chem.BondType.SINGLE:
            return False, "Unsaturated carbon-carbon bond found"
    
    # Get the carbonyl carbon and start chain traversal
    matches = mol.GetSubstructMatches(carboxylic_acid)
    carbonyl_idx = matches[0][0]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    
    # Find adjacent carbon in the main chain
    chain_next = None
    for neighbor in carbonyl_atom.GetNeighbors():
        if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() not in (matches[0][1], matches[0][2]):
            chain_next = neighbor
            break
    if not chain_next:
        return False, "No chain attached to carboxylic acid"
    
    # Traverse main chain to check for branches
    chain_carbons = {carbonyl_idx}
    current = chain_next
    prev_idx = carbonyl_idx
    
    while True:
        if current.GetIdx() in chain_carbons:  # Prevent infinite loops
            break
        chain_carbons.add(current.GetIdx())
        
        # Get next carbons in chain (exclude previous)
        next_carbons = []
        for neighbor in current.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() != prev_idx:
                next_carbons.append(neighbor)
        
        # Branch check
        if len(next_carbons) > 1:
            return False, "Branch detected in carbon chain"
        if not next_carbons:
            break  # End of chain
        
        prev_idx = current.GetIdx()
        current = next_carbons[0]
    
    # Verify all carbons are in the main chain
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetIdx() not in chain_carbons:
            return False, "Extra carbon detected (branch or substituent)"
    
    return True, "Straight-chain saturated fatty acid with carboxylic acid group"