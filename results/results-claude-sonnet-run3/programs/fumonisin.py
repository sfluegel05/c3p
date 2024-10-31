from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fumonisin(smiles: str):
    """
    Determines if a molecule is a fumonisin.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a fumonisin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    
    # Check for key structural features
    # - Long carbon chain (at least 18-20 carbons)
    # - Multiple carboxylic acid groups
    # - Amino group
    # - Multiple hydroxyl groups
    # - Ester linkages
    
    # Count key functional groups
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)O')
    ester_pattern = Chem.MolFromSmarts('C(=O)OC')
    amino_pattern = Chem.MolFromSmarts('[NH2,NH]')
    hydroxyl_pattern = Chem.MolFromSmarts('[OH]')
    
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    num_carboxylic = len(mol.GetSubstructMatches(carboxylic_acid_pattern))
    num_esters = len(mol.GetSubstructMatches(ester_pattern))
    num_amino = len(mol.GetSubstructMatches(amino_pattern))
    num_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))

    # Check structural requirements
    if num_carbons < 18:
        return False, "Carbon chain too short for fumonisin"
        
    if num_carboxylic < 2:
        return False, "Insufficient carboxylic acid groups"
        
    if num_esters < 2:
        return False, "Insufficient ester linkages"
        
    if num_amino < 1:
        return False, "Missing amino group"
        
    if num_hydroxyls < 2:
        return False, "Insufficient hydroxyl groups"

    # If all requirements are met, likely a fumonisin
    return True, f"Matches fumonisin structure: {num_carbons} carbons, {num_carboxylic} COOH groups, {num_esters} esters, {num_amino} NH2 groups, {num_hydroxyls} OH groups"
# Pr=1.0
# Recall=1.0