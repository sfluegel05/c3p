"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a broader substructure pattern for the steroid core including stereochemistry
    steroid_patterns = [
        Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CC(C4=C3CC[C@@H]4)C'), # A possible general steroid pattern
        Chem.MolFromSmarts('C1CC[C@H]2[C@@H]([C@H]1)CC[C@@H]1[C@H]2CC[C@@H]1C') # Including stereochemistry
    ]
    
    if not any(mol.HasSubstructMatch(steroid_pattern) for steroid_pattern in steroid_patterns):
        return False, "No steroid core found"

    # Check for glycosidic linkages in a broader sense
    glycoside_patterns = [
        Chem.MolFromSmarts('OC1C(O)C(O)C(O)C(O)C1O'),
        Chem.MolFromSmarts('OC1OCC(O)C(O)C1O'),  # Alternate sugar ring patterns
    ]
    
    if not any(mol.HasSubstructMatch(glycoside_pattern) for glycoside_pattern in glycoside_patterns):
        return False, "No glycosidic linkage found"

    # Check for the presence of a lactone ring or similar ester structures
    lactone_patterns = [
        Chem.MolFromSmarts('O=C1COCC1'),  # Common lactone structure
        Chem.MolFromSmarts('O=C1CCC(O1)') # Broader cyclic ester patterns
    ]
    
    if not any(mol.HasSubstructMatch(lactone_pattern) for lactone_pattern in lactone_patterns):
        return False, "No lactone or related ester ring found"
    
    # Additional characterizations like hydroxyl groups
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetTotalDegree() == 1)
    if hydroxyl_count < 3:
        return False, f"Insufficient hydroxyl groups in the molecule: {hydroxyl_count} found"
    
    return True, "The structure contains a steroid core with glycosidic linkage and a lactone ring characteristic of cardiac glycosides"

# Example usage
# result = is_cardiac_glycoside("SMILES_STRING_HERE")
# print(result)