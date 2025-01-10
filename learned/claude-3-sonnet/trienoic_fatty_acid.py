"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: trienoic fatty acid
Definition: Any polyunsaturated fatty acid that contains three double bonds.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count number of double bonds
    double_bond_pattern = Chem.MolFromSmarts('C=C')
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    
    # Also count cumulated double bonds (e.g., allenes C=C=C)
    allene_pattern = Chem.MolFromSmarts('C=C=C')
    allene_bonds = len(mol.GetSubstructMatches(allene_pattern))
    
    # Each allene contains two double bonds, so adjust total count
    total_double_bonds = double_bonds - allene_bonds
    
    if total_double_bonds != 3:
        return False, f"Found {total_double_bonds} double bonds, need exactly 3"
    
    # Check carbon chain length (typical fatty acids have >8 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 8:
        return False, f"Carbon chain too short ({carbon_count} carbons, need at least 8)"
    
    # Check for reasonable molecular weight (typical fatty acids are 200-400 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 400:
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for fatty acids (200-400 Da)"
    
    # Check for continuous carbon chain
    alkyl_chain_pattern = Chem.MolFromSmarts('CCCC')
    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "No continuous carbon chain found"
    
    # Additional check to exclude complex molecules that happen to have 3 double bonds
    # Count number of rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 2:  # Allow up to 2 rings for epoxy fatty acids
        return False, f"Too many rings ({ring_count}) for a fatty acid"
    
    # Count number of oxygen atoms (should be reasonable for a fatty acid)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count > 6:  # Allow for hydroxy and epoxy groups
        return False, f"Too many oxygen atoms ({oxygen_count}) for a fatty acid"
        
    return True, "Contains carboxylic acid group and exactly three double bonds in a fatty acid structure"