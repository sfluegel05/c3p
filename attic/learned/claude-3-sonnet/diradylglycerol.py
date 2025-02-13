"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: diradylglycerol
A lipid with glycerol bearing two substituent groups at any two of three possible positions
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern (C-C-C)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for the specific pattern of glycerol with attachments
    # [O,C] means oxygen or carbon can be attached to the glycerol carbons
    glycerol_with_attachments = Chem.MolFromSmarts("[CH2X4]([O,C])[CHX4]([O,C])[CH2X4]([O,C])")
    if not mol.HasSubstructMatch(glycerol_with_attachments):
        return False, "Glycerol backbone lacks proper substitution pattern"
    
    # Count free hydroxyl groups on the glycerol backbone
    # This pattern looks for -OH groups connected to the glycerol carbons
    hydroxyl_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-[OX2H]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_matches == 0:
        return False, "No free hydroxyl groups found on glycerol backbone"
    
    if hydroxyl_matches > 1:
        return False, "Too many free hydroxyl groups for diradylglycerol"
        
    # Look for ester linkages specifically attached to glycerol backbone
    ester_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-[OX2][CX3](=[OX1])[#6]")
    ester_count = len(mol.GetSubstructMatches(ester_pattern))
    
    # Look for ether linkages specifically attached to glycerol backbone
    ether_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-[OX2]-[CX4]")
    ether_count = len(mol.GetSubstructMatches(ether_pattern))
    
    total_substituents = ester_count + ether_count
    
    if total_substituents != 2:
        return False, f"Found {total_substituents} substituents attached to glycerol, need exactly 2"
    
    # Verify presence of long carbon chains
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = len(mol.GetSubstructMatches(chain_pattern))
    if chain_matches < 2:
        return False, "Missing long carbon chains"
    
    # Count carbons to verify chains are long enough
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Carbon chains too short for diradylglycerol"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for diradylglycerol"
        
    return True, "Contains glycerol backbone with two substituent groups and one free hydroxyl"