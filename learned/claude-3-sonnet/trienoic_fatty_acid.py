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
    if not mol.HasSubstructMatches(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carbon chain length from carboxylic acid
    long_chain_pattern = Chem.MolFromSmarts('C(=O)[O][H].{8,}')
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Carbon chain too short for a fatty acid"
    
    # Count rings but exclude small rings (epoxides)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    epoxide_pattern = Chem.MolFromSmarts('C1OC1')
    epoxide_count = len(mol.GetSubstructMatches(epoxide_pattern))
    non_epoxide_rings = ring_count - epoxide_count
    
    if non_epoxide_rings > 0:
        return False, f"Contains {non_epoxide_rings} non-epoxide rings"
    
    # Count double bonds more accurately
    double_bond_pattern = Chem.MolFromSmarts('C=C')
    allene_pattern = Chem.MolFromSmarts('C=C=C')
    cumulated_pattern = Chem.MolFromSmarts('C=C=C=C')
    
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    allenes = len(mol.GetSubstructMatches(allene_pattern))
    cumulated = len(mol.GetSubstructMatches(cumulated_pattern))
    
    # Adjust count for cumulated systems
    total_double_bonds = double_bonds - allenes + (allenes - cumulated)
    
    if total_double_bonds != 3:
        return False, f"Found {total_double_bonds} double bonds, need exactly 3"
    
    # Check for reasonable fatty acid characteristics
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:  # Most natural fatty acids have at least 12 carbons
        return False, f"Carbon chain too short ({carbon_count} carbons)"
    
    # Check for continuous aliphatic chain
    aliphatic_chain = Chem.MolFromSmarts('CCCCC')  # At least 5 consecutive carbons
    if not mol.HasSubstructMatch(aliphatic_chain):
        return False, "No significant aliphatic chain found"
    
    # Check for aromatic rings
    aromatic_pattern = Chem.MolFromSmarts('a')
    if mol.HasSubstructMatch(aromatic_pattern):
        return False, "Contains aromatic rings"
    
    # Allow hydroxyl, epoxy, and hydroperoxy groups common in fatty acids
    if ring_count > 2:  # Allow maximum 2 epoxide groups
        return False, f"Too many rings ({ring_count})"
        
    return True, "Contains carboxylic acid group and exactly three double bonds in a fatty acid structure"