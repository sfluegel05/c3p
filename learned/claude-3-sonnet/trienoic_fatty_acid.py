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
        # Check if it's part of a larger lipid (ester)
        ester_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2]')
        if not mol.HasSubstructMatch(ester_pattern):
            return False, "No carboxylic acid or ester group found"
    
    # Count carbon chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:  # Most natural fatty acids have at least 12 carbons
        return False, f"Carbon chain too short ({carbon_count} carbons)"
    
    # Look for main fatty acid chain pattern
    fatty_chain = Chem.MolFromSmarts('CCCCC')  # At least 5 consecutive carbons
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "No significant aliphatic chain found"
    
    # Count double bonds more carefully
    double_bond_patterns = [
        ('C=C', 'standard'),
        ('C/C=C/C', 'trans'),
        ('C/C=C\\C', 'cis')
    ]
    
    total_double_bonds = 0
    for pattern, _ in double_bond_patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(pattern_mol)
        total_double_bonds += len(matches)
    
    # Normalize double bond count (avoid double counting)
    total_double_bonds = total_double_bonds // 2
    
    if total_double_bonds != 3:
        return False, f"Found {total_double_bonds} double bonds, need exactly 3"
    
    # Check for aromatic rings (not allowed in fatty acids)
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"
    
    # Allow common substituents but limit their number
    allowed_groups = [
        ('[OH]', 'hydroxyl'),
        ('[OOH]', 'hydroperoxy'),
        ('C1OC1', 'epoxy'),
        ('C(=O)C', 'keto')
    ]
    
    substituent_count = 0
    for pattern, name in allowed_groups:
        pattern_mol = Chem.MolFromSmarts(pattern)
        matches = len(mol.GetSubstructMatches(pattern_mol))
        substituent_count += matches
    
    if substituent_count > 3:
        return False, f"Too many substituents ({substituent_count}) for a typical fatty acid"
    
    # Check for reasonable oxygen count
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count > 8:
        return False, f"Too many oxygen atoms ({oxygen_count}) for a fatty acid"
    
    # Verify chain connectivity (avoid highly branched structures)
    chain_pattern = Chem.MolFromSmarts('CCCCCCCC')  # Look for 8 consecutive carbons
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No continuous carbon chain of sufficient length"
        
    return True, "Contains a fatty acid structure with exactly three double bonds"