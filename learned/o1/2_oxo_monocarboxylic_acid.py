"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: CHEBI:37727 2-oxo monocarboxylic acid
"""
from rdkit import Chem

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
    
    # Define SMARTS patterns
    # Pattern for monocarboxylic acid group (protonated or deprotonated)
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H0-,OX2H1]')
    # Pattern for carbonyl group attached to carbon (ketone)
    oxo_group_pattern = Chem.MolFromSmarts('[CX3]=O')
    
    # Find all carboxylic acid groups
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) == 0:
        return False, "No carboxylic acid groups found"
    elif len(carboxylic_acid_matches) > 1:
        return False, f"More than one carboxylic acid group found ({len(carboxylic_acid_matches)} groups)"
    
    # There is exactly one carboxylic acid group
    acid_carbon_idx = carboxylic_acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Find alpha carbons (carbons adjacent to the carboxylic acid carbon)
    alpha_carbons = [neighbor for neighbor in acid_carbon.GetNeighbors() if neighbor.GetAtomicNum() == 6]
    
    if not alpha_carbons:
        return False, "No alpha carbon adjacent to carboxylic acid carbon"
    
    # Check if any alpha carbon has an oxo substituent
    found_2_oxo = False
    for alpha_carbon in alpha_carbons:
        alpha_carbon_idx = alpha_carbon.GetIdx()
        # Check for oxo group attached to alpha carbon
        for neighbor in alpha_carbon.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(alpha_carbon_idx, neighbor.GetIdx())
            if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                found_2_oxo = True
                break
        if found_2_oxo:
            break
    
    if not found_2_oxo:
        return False, "No 2-oxo substituent found on alpha carbon"
    
    return True, "Molecule is a 2-oxo monocarboxylic acid"