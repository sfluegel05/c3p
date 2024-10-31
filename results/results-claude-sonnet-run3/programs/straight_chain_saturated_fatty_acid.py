from rdkit import Chem
from rdkit.Chem import AllChem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight chain saturated fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a straight chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
        
    # Count number of carboxylic acid groups
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(matches) > 1:
        return False, "Multiple carboxylic acid groups found"

    # Get the carboxyl carbon atom
    carboxyl_carbon = mol.GetAtomWithIdx(matches[0][0])
    
    # Check for branching
    for atom in mol.GetAtoms():
        # Skip carboxyl group atoms
        if atom.GetIdx() in matches[0]:
            continue
        # Allow hydroxy substituents
        if atom.GetSymbol() == 'O' and atom.GetTotalDegree() == 2:
            continue
        # Check for branching - atoms should have at most 2 carbon neighbors
        carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() in ['C','[2H]']]
        if len(carbon_neighbors) > 2:
            return False, "Branched carbon chain found"
            
    # Check for unsaturation by looking for double/triple bonds
    # Skip the carboxyl C=O bond
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() > 1:
            # Skip the carboxyl C=O
            if bond.GetBeginAtomIdx() in matches[0] and bond.GetEndAtomIdx() in matches[0]:
                continue
            return False, "Contains unsaturated bonds"
        
    # Count carbons in main chain
    carbon_count = len([a for a in mol.GetAtoms() if a.GetSymbol() == 'C'])
    
    return True, f"Straight chain saturated fatty acid with {carbon_count} carbons"
# Pr=1.0
# Recall=0.975609756097561