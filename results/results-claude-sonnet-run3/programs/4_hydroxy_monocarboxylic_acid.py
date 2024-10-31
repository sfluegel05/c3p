from rdkit import Chem
from rdkit.Chem import AllChem

def is_4_hydroxy_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 4-hydroxy monocarboxylic acid (has a hydroxy group gamma to the carboxy group).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 4-hydroxy monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Find carboxylic acid groups
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"
    
    if len(carboxyl_matches) > 1:
        return False, "More than one carboxylic acid group found"

    # Find hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts('[CH,CH2,CH3][OH]')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if not hydroxyl_matches:
        return False, "No hydroxyl group found"

    # Get the carbon atom of the carboxyl group
    carboxyl_carbon = carboxyl_matches[0][0]
    
    # For each hydroxyl group, check if it's gamma (4 carbons away) to the carboxyl group
    for hydroxyl_match in hydroxyl_matches:
        hydroxyl_carbon = hydroxyl_match[0]
        
        # Find shortest path between hydroxyl carbon and carboxyl carbon
        path = Chem.GetShortestPath(mol, hydroxyl_carbon, carboxyl_carbon)
        
        if path is None:
            continue
            
        # Path length should be 4 for gamma position (C-C-C-C(=O)OH)
        if len(path) == 4:
            # Verify the path consists of carbon chain
            atoms_in_path = [mol.GetAtomWithIdx(i).GetSymbol() for i in path]
            if all(atom == 'C' for atom in atoms_in_path):
                return True, "Found hydroxyl group gamma to carboxylic acid"

    return False, "No hydroxyl group found in gamma position to carboxylic acid"
# Pr=1.0
# Recall=0.9166666666666666