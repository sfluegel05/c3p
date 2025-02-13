"""
Classifies: CHEBI:35366 fatty acid
"""
"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carboxylic acid groups - should only have one
    carboxyl_matches = len(mol.GetSubstructMatches(carboxyl_pattern))
    if carboxyl_matches > 1:
        return False, f"Found {carboxyl_matches} carboxylic acid groups, should have exactly 1"
    
    # Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 2:  # Minimum 2 carbons (one for COOH)
        return False, f"Too few carbons ({carbon_count}) for a fatty acid"
    if carbon_count > 35:  # Allow some margin above typical max of 28
        return False, f"Too many carbons ({carbon_count}) for a fatty acid"
        
    # Check if the molecule is entirely aliphatic (excluding the carboxyl group)
    aromatic_pattern = Chem.MolFromSmarts("a")
    if mol.HasSubstructMatch(aromatic_pattern):
        # Some exceptions exist where a small aromatic group might be present
        # Count aromatic atoms
        aromatic_atoms = len(mol.GetSubstructMatches(aromatic_pattern))
        if aromatic_atoms > 6:  # Allow for one benzene ring
            return False, "Contains too many aromatic atoms"
    
    # Look for a carbon chain of at least 3 carbons (plus the carboxyl carbon)
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No sufficient carbon chain found"
    
    # Additional checks for reasonable molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 60:  # Less than acetic acid
        return False, "Molecular weight too low for fatty acid"
    if mol_wt > 500:  # Most fatty acids are under 500 Da
        # Some exceptions for heavily substituted fatty acids
        if mol_wt > 1000:
            return False, "Molecular weight too high for fatty acid"
    
    # Check the ratio of carbons to other atoms
    # Fatty acids typically have a high C:X ratio
    total_atoms = sum(1 for atom in mol.GetAtoms())
    if carbon_count / total_atoms < 0.4:  # Less than 40% carbon
        return False, "Too few carbons relative to other atoms"
    
    # If we've made it here, classify based on features
    features = []
    
    # Check for unsaturation
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if mol.HasSubstructMatch(double_bond_pattern):
        features.append("unsaturated")
    if mol.HasSubstructMatch(triple_bond_pattern):
        features.append("contains triple bonds")
    
    # Check for substitution
    if any(atom.GetAtomicNum() not in [1,6,8] for atom in mol.GetAtoms()):
        features.append("substituted")
    
    # Check for hydroxy groups (besides COOH)
    hydroxy_pattern = Chem.MolFromSmarts("[CX4][OX2H1]")
    if mol.HasSubstructMatch(hydroxy_pattern):
        features.append("hydroxylated")
    
    feature_str = " and ".join(features) if features else "saturated"
    return True, f"Aliphatic carboxylic acid with {carbon_count} carbons, {feature_str}"