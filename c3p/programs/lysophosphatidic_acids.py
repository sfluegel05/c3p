"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: CHEBI:75856 lysophosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phosphate group (-OP(=O)(O)O)
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[O])([OX2H,OX1-])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Look for glycerol backbone pattern (C-C-C with oxygens)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for single ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"
    
    # Look for fatty acid chain (long carbon chain attached to ester)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No fatty acid chain found"
    
    # Count rotatable bonds to verify chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
        return False, "Chain too short to be a fatty acid"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 5:
        return False, "Too few carbons for lysophosphatidic acid"
    if o_count < 6:
        return False, "Must have at least 6 oxygens"
    if p_count != 1:
        return False, "Must have exactly one phosphorus"
        
    # Look for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 2:
        return False, "Must have at least 2 hydroxyl groups"
    
    return True, "Contains glycerol backbone with one fatty acid chain and phosphate group"