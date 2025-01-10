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
    
    # Count carbon chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:  # Most natural fatty acids have at least 12 carbons
        return False, f"Carbon chain too short ({carbon_count} carbons)"
    
    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts('C=C')
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    
    # Check for triple bonds (should not be present)
    triple_bond_pattern = Chem.MolFromSmarts('C#C')
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Contains triple bonds"
    
    if double_bonds != 3:
        return False, f"Found {double_bonds} double bonds, need exactly 3"
    
    # Check for aromatic rings
    aromatic_pattern = Chem.MolFromSmarts('a')
    if mol.HasSubstructMatch(aromatic_pattern):
        return False, "Contains aromatic rings"
    
    # Count rings but exclude small rings (epoxides)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    epoxide_pattern = Chem.MolFromSmarts('C1OC1')
    epoxide_matches = len(mol.GetSubstructMatches(epoxide_pattern))
    non_epoxide_rings = ring_count - epoxide_matches
    
    if non_epoxide_rings > 0:
        return False, f"Contains {non_epoxide_rings} non-epoxide rings"
        
    # Allow common substituents found in fatty acids
    allowed_groups = [
        ('[OH]', 'hydroxyl'),
        ('[OOH]', 'hydroperoxy'),
        ('C1OC1', 'epoxy')
    ]
    
    # Check if molecule has a reasonable fatty acid backbone
    aliphatic_chain = Chem.MolFromSmarts('CCCC')  # At least 4 consecutive carbons
    if not mol.HasSubstructMatch(aliphatic_chain):
        return False, "No significant aliphatic chain found"
    
    # Verify oxygen count is reasonable (carboxylic acid = 2, plus allowed substituents)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count > 8:  # Allow up to 3 additional oxygen-containing groups
        return False, f"Too many oxygen atoms ({oxygen_count}) for a fatty acid"
        
    return True, "Contains carboxylic acid group and exactly three double bonds in a fatty acid structure"