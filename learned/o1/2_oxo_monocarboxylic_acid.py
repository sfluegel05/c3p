"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: CHEBI:37727 2-oxo monocarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid is any monocarboxylic acid having a 2-oxo substituent.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for monocarboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H]')
    acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    # Check for exactly one monocarboxylic acid group
    if len(acid_matches) != 1:
        return False, f"Found {len(acid_matches)} carboxylic acid groups, need exactly 1"
    
    # Get the index of the carboxylic acid carbon
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Get neighboring atoms of the carboxylic acid carbon (excluding oxygen atoms in the acid group)
    acid_group_indices = set(acid_matches[0])
    alpha_carbons = []
    for neighbor in acid_carbon.GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx not in acid_group_indices and neighbor.GetAtomicNum() == 6:
            alpha_carbons.append(neighbor)
    
    # Check if there is exactly one alpha carbon
    if len(alpha_carbons) != 1:
        return False, f"Found {len(alpha_carbons)} alpha carbons, need exactly 1"
    
    alpha_carbon = alpha_carbons[0]
    
    # Check if the alpha carbon has a ketone (=O) group
    has_alpha_oxo = False
    for bond in alpha_carbon.GetBonds():
        neighbor = bond.GetOtherAtom(alpha_carbon)
        if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            has_alpha_oxo = True
            break
    
    if not has_alpha_oxo:
        return False, "Alpha carbon does not have a ketone (=O) group"
    
    return True, "Molecule is a 2-oxo monocarboxylic acid"