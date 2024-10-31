from rdkit import Chem
from rdkit.Chem import AllChem

def is_fluoroalkanoic_acid(smiles: str):
    """
    Determines if a molecule is a fluoroalkanoic acid (perfluorinated derivative of alkanoic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fluoroalkanoic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
        
    # Get the carboxylic carbon atom
    matches = mol.GetSubstructMatches(carboxylic_pattern)
    if not matches:
        return False, "Could not identify carboxylic carbon"
    carboxylic_carbon_idx = matches[0][0]
    
    # Find all carbons in the molecule
    carbon_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            carbon_atoms.append(atom)
            
    # Check if all carbons (except carboxylic) are fully fluorinated
    for atom in carbon_atoms:
        if atom.GetIdx() == carboxylic_carbon_idx:
            continue
            
        # Count fluorine neighbors
        neighbors = atom.GetNeighbors()
        f_count = sum(1 for n in neighbors if n.GetSymbol() == 'F')
        c_count = sum(1 for n in neighbors if n.GetSymbol() == 'C')
        
        # Each non-carboxylic carbon should have exactly 2 fluorines and 2 carbon neighbors
        # or 3 fluorines and 1 carbon neighbor (for internal carbons)
        # or all fluorines for terminal carbon
        if not ((f_count == 2 and c_count == 2) or 
                (f_count == 3 and c_count == 1) or
                (f_count == atom.GetTotalValence() and c_count == 0)):
            return False, f"Carbon at position {atom.GetIdx()} is not properly fluorinated"
            
    chain_length = len(carbon_atoms)
    if chain_length < 2:
        return False, "Chain too short for a fluoroalkanoic acid"
        
    return True, f"Perfluorinated alkanoic acid with chain length {chain_length}"
# Pr=1.0
# Recall=1.0