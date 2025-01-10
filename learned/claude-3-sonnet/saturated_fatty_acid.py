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
    
    # Count carboxylic acid groups (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxyl_matches = len(mol.GetSubstructMatches(carboxyl_pattern))
    if carboxyl_matches == 0:
        return False, "No carboxylic acid group found"
    if carboxyl_matches > 1:
        return False, "Multiple carboxylic acid groups found"
    
    # Check for carbon-carbon multiple bonds
    cc_double_bond = Chem.MolFromSmarts("[C]=[C;!$(C(=O)[OH])]")
    if mol.HasSubstructMatch(cc_double_bond):
        return False, "Contains carbon-carbon double bonds"
    
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Contains carbon-carbon triple bonds"
    
    # Check for aromatic systems
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic rings"
    
    # Count carbons and check atom types
    c_count = 0
    o_count = 0
    other_atoms = False
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:  # Carbon
            c_count += 1
        elif atomic_num == 8:  # Oxygen
            o_count += 1
        elif atomic_num == 1:  # Hydrogen (including deuterium)
            continue
        else:
            other_atoms = True
    
    if other_atoms:
        return False, "Contains atoms other than C, H, O"
    
    # Check oxygen count - fatty acids should have 2 oxygens (from COOH)
    # Allow up to 3 oxygens to account for hydroxy fatty acids
    if o_count > 3:
        return False, "Too many oxygen atoms for a fatty acid"
    
    # Check for cyclic structures (except small rings like cyclopropane)
    ring_info = mol.GetRingInfo()
    for ring_size in ring_info.RingSizes():
        if ring_size > 3:
            return False, "Contains large ring structures"
    
    # Check for excessive branching
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            if len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]) > 3:
                return False, "Contains quaternary carbon centers"
    
    # Verify chain length and structure
    if c_count < 2:
        return False, "Carbon chain too short for fatty acid"
        
    # Look for significant aliphatic chain
    # More lenient pattern that allows branching
    aliphatic_chain = Chem.MolFromSmarts("CC")
    if not mol.HasSubstructMatch(aliphatic_chain):
        return False, "No aliphatic chain found"
    
    # Additional check for ketone groups (which shouldn't be present)
    ketone_pattern = Chem.MolFromSmarts("[C!$(C(=O)O)]=O")
    if mol.HasSubstructMatch(ketone_pattern):
        return False, "Contains ketone groups"
    
    # Check for ester groups
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if mol.HasSubstructMatch(ester_pattern):
        return False, "Contains ester groups"
        
    return True, "Saturated fatty acid with aliphatic chain and no carbon-carbon multiple bonds"