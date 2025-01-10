"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: saturated fatty acids
Definition: Any fatty acid containing no carbon to carbon multiple bonds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Check for carboxylic acid group (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for carbon-carbon double bonds
    # Specifically exclude the carboxylic acid carbon
    cc_double_bond = Chem.MolFromSmarts("[C]=[C;!$(C(=O)[OH])]")
    if mol.HasSubstructMatch(cc_double_bond):
        return False, "Contains carbon-carbon double bonds"
    
    # Check for carbon-carbon triple bonds
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Contains carbon-carbon triple bonds"
    
    # Check for aromatic systems
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"
    
    # Count carbons in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 2:
        return False, "Carbon chain too short for fatty acid"
    
    # Check for aliphatic chain
    aliphatic_chain = Chem.MolFromSmarts("CCCC")  # At least 4 carbons in chain
    if not mol.HasSubstructMatch(aliphatic_chain):
        return False, "No significant aliphatic chain found"
    
    # Check atom types (allowing for deuterium)
    allowed_atoms = {1, 6, 8}  # H, C, O
    atom_nums = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if not atom_nums.issubset(allowed_atoms):
        return False, "Contains atoms other than C, H, O"
        
    return True, "Saturated fatty acid with aliphatic chain and no carbon-carbon multiple bonds"