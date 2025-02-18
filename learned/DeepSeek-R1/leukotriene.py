"""
Classifies: CHEBI:25029 leukotriene
"""
"""
Classifies: CHEBI:24400 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    Leukotrienes are C20 icosanoids with four double bonds, three of which are conjugated.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for exactly four double bonds
    double_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE]
    if len(double_bonds) != 4:
        return False, f"Found {len(double_bonds)} double bonds (need exactly 4)"
    
    # Check for conjugated triene (three consecutive conjugated double bonds)
    triene_pattern = Chem.MolFromSmarts('[*]=[*]-[*]=[*]-[*]=[*]')
    if not mol.HasSubstructMatch(triene_pattern):
        return False, "No conjugated triene found"
    
    # Check carbon count (C20 backbone)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, f"Only {carbon_count} carbons (minimum 20 required)"
    
    # Optional: Check for carboxylic acid group (common but not universal)
    # carboxylic_acid = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    # if not mol.HasSubstructMatch(carboxylic_acid):
    #     return False, "No carboxylic acid group"
    
    return True, "C20 with four double bonds (three conjugated)"