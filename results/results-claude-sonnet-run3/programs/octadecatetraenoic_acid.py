from rdkit import Chem
from rdkit.Chem import AllChem

def is_octadecatetraenoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecatetraenoic acid (18-carbon chain with 4 double bonds and a carboxylic acid group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecatetraenoic acid, False otherwise
        str: Reason for classification
    """
    # Handle invalid SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count total carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count != 18:
        return False, f"Carbon count is {carbon_count}, should be 18"

    # Count all double bonds
    double_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Don't count the carboxylic acid double bond
            if not (bond.GetBeginAtom().GetSymbol() == 'O' or bond.GetEndAtom().GetSymbol() == 'O'):
                double_bonds.append(bond)

    if len(double_bonds) != 4:
        return False, f"Carbon-carbon double bond count is {len(double_bonds)}, should be 4"

    # Check if it's a straight chain
    # Get the carboxylic acid carbon
    carb_pattern = Chem.MolFromSmarts('C(=O)O')
    matches = mol.GetSubstructMatches(carb_pattern)
    if not matches:
        return False, "Could not identify carboxylic acid carbon"
    
    carb_carbon = matches[0][0]
    
    # Traverse the chain
    visited = set()
    current = carb_carbon
    chain = []
    
    while True:
        visited.add(current)
        chain.append(current)
        
        # Get non-oxygen neighbors we haven't visited
        next_atoms = [n.GetIdx() for n in mol.GetAtomWithIdx(current).GetNeighbors() 
                     if n.GetSymbol() != 'O' and n.GetIdx() not in visited]
        
        if not next_atoms:
            break
            
        # Check for branching
        if len(next_atoms) > 1:
            carbon_neighbors = [n for n in next_atoms 
                              if mol.GetAtomWithIdx(n).GetSymbol() == 'C']
            if len(carbon_neighbors) > 1:
                return False, "Molecule contains branching - not a straight chain"
        
        current = next_atoms[0]

    # Verify chain length
    carbon_chain = [i for i in chain if mol.GetAtomWithIdx(i).GetSymbol() == 'C']
    if len(carbon_chain) != 18:
        return False, f"Main carbon chain length is {len(carbon_chain)}, should be 18"

    return True, "Valid octadecatetraenoic acid with 18 carbons, 4 double bonds in straight chain"
# Pr=1.0
# Recall=1.0