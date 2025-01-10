"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: olefinic fatty acid
Definition: Any fatty acid containing at least one C=C double bond
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for carbon-carbon double bond
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No carbon-carbon double bond found"
    
    # Count carbons to ensure it's a fatty acid (typically >4 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4:
        return False, "Carbon chain too short to be a fatty acid"
    
    # Additional checks to avoid false positives
    
    # Check for continuous carbon chain
    aliphatic_chain = Chem.MolFromSmarts("[C]~[C]~[C]~[C]")
    if not mol.HasSubstructMatch(aliphatic_chain):
        return False, "No continuous carbon chain found"
    
    # Count rotatable bonds to verify chain flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Too rigid to be a fatty acid"
    
    # Get the number of double bonds
    double_bond_matches = len(mol.GetSubstructMatches(double_bond_pattern))
    
    # Calculate degree of unsaturation
    rings = rdMolDescriptors.CalcNumRings(mol)
    
    reason = f"Contains carboxylic acid group and {double_bond_matches} C=C double bond"
    if double_bond_matches > 1:
        reason += "s"
    if rings > 0:
        reason += f" ({rings} ring structure{'s' if rings > 1 else ''} present)"
    
    return True, reason