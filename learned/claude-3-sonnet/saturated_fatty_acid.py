"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: saturated fatty acids
Definition: Any fatty acid containing no carbon to carbon multiple bonds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (-COOH) or ester group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C,H]")
    
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    
    if not (has_carboxyl or has_ester):
        return False, "No carboxylic acid or ester group found"
    
    # Check for absence of carbon-carbon multiple bonds
    multiple_bond_pattern = Chem.MolFromSmarts("[C]=[C,N,O;!R]")
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    
    if mol.HasSubstructMatch(multiple_bond_pattern):
        return False, "Contains carbon-carbon double bonds"
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Contains carbon-carbon triple bonds"
    
    # Check for aromatic systems
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"
    
    # Count rings - fatty acids should have minimal rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 2:  # Allow up to 2 rings for cyclopropane fatty acids
        return False, "Too many ring structures"
        
    # Check for sugar-like patterns (multiple OH groups in close proximity)
    polyol_pattern = Chem.MolFromSmarts("[OH]-[CH]-[CH]-[CH]-[OH]")
    if mol.HasSubstructMatch(polyol_pattern):
        return False, "Contains sugar-like hydroxyl pattern"
    
    # Count carbons and check chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 2:
        return False, "Carbon chain too short for fatty acid"
    
    # Check for long aliphatic chain
    aliphatic_chain = Chem.MolFromSmarts("CCCC")  # At least 4 carbons in chain
    if not mol.HasSubstructMatch(aliphatic_chain):
        return False, "No significant aliphatic chain found"
    
    # Check atom types (allowing for deuterium)
    allowed_atoms = {1, 6, 8}  # H, C, O
    atom_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if not atom_nums.issubset(allowed_atoms):
        return False, "Contains atoms other than C, H, O"
    
    # Count ketone groups (excluding the carboxylic acid carbonyl)
    ketone_pattern = Chem.MolFromSmarts("[C!$(C(=O)O)]=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if len(ketone_matches) > 0:
        return False, "Contains additional ketone groups"
    
    # Calculate fraction of sp3 carbons (should be high for saturated fatty acids)
    sp3_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4]")))
    if sp3_carbons / c_count < 0.7:  # At least 70% of carbons should be sp3
        return False, "Insufficient sp3 carbons for saturated fatty acid"
    
    return True, "Saturated fatty acid or derivative with aliphatic chain"