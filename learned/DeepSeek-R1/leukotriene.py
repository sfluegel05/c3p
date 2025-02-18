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
    Leukotrienes are C20 icosanoids with four double bonds (three conjugated in a triene structure).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check C20 backbone (at least 20 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, f"Only {carbon_count} carbons (need ≥20)"
    
    # Get all double bonds (including those in functional groups)
    double_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE]
    if len(double_bonds) < 4:
        return False, f"Found {len(double_bonds)} double bonds (need ≥4)"
    
    # Check for conjugated triene (three consecutive conjugated double bonds)
    # Updated pattern to match any three consecutive conjugated double bonds
    triene_pattern = Chem.MolFromSmarts('C=CC=CC=C')
    if not mol.HasSubstructMatch(triene_pattern):
        return False, "No conjugated triene found"
    
    # Check oxygen content (leukotrienes typically have multiple oxygen atoms)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3:
        return False, f"Only {oxygen_count} oxygens (common in leukotrienes)"
    
    # Check for characteristic carboxylic acid/derivative (ester/amide common)
    carbonyl_group = Chem.MolFromSmarts('[CX3]=[OX1]')
    if not mol.HasSubstructMatch(carbonyl_group):
        return False, "No carbonyl group found"
    
    return True, "C20 with ≥4 double bonds (three consecutive conjugated triene), oxygen content, and carbonyl group"