"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies chemical entities of the class CHEBI:17514 quinic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid based on its SMILES string.
    A quinic acid is defined as a cyclitol carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyclohexane core
    cyclohexane_pattern = Chem.MolFromSmarts("[C&R1&R2&R3&R4&R5&R6]1[C&R2&R3&R4&R5&R6][C&R3&R4&R5&R6][C&R4&R5&R6][C&R5&R6][C&R6]1")
    if not mol.HasSubstructMatch(cyclohexane_pattern):
        return False, "No cyclohexane core found"
    
    # Check for at least 3 hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, "Fewer than 3 hydroxyl groups found"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OX2H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, "Did not find exactly 1 carboxylic acid group"
    
    # Check for specific stereochemistry of quinic acid
    quinic_acid_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@@H]([C@@H](O)O)O)O)O)O)C(=O)O")
    if not mol.HasSubstructMatch(quinic_acid_pattern):
        return False, "Stereochemistry does not match quinic acid"
    
    # Check for common substituents
    allowed_substituents = ["[OC]", "[C&R1&R2&R3&R4&R5&R6]", "[C&R1&R2&R3&R4&R5&R6][C&R1&R2&R3&R4&R5&R6]", "[C](=O)[OX2H1]", "[C](=O)[OX2][C&R1&R2&R3&R4&R5&R6]"]
    substituents = [Chem.MolToSmarts(frag) for frag in Chem.GetMolFrags(mol, aromatics=Chem.AromaticityMatchers.AromaticAny())]
    for sub in substituents:
        if sub not in allowed_substituents:
            return False, f"Unexpected substituent found: {sub}"
    
    # Check molecular weight (quinic acids typically 192-500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 192 or mol_wt > 500:
        return False, "Molecular weight outside typical range for quinic acids"
    
    return True, "Molecule matches the structure and properties of a quinic acid"