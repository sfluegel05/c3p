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
    
    # Look for glycerol backbone pattern (C-C-C with oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Look for ether groups (-O-C-)
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    # Count total substituents (esters + ethers)
    # Note: Need to filter out the glycerol backbone oxygens from ether count
    total_substituents = len(ester_matches) + len(ether_matches)
    
    # Count hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # For diradylglycerols, we expect:
    # - Exactly 2 substituents (ester and/or ether)
    # - At least 1 hydroxyl group
    # - Long carbon chains attached to the substituents
    
    if len(ester_matches) + len(ether_matches) != 2:
        return False, f"Found {total_substituents} substituents, need exactly 2"
    
    if len(hydroxyl_matches) < 1:
        return False, "No free hydroxyl group found"
    
    # Check for long carbon chains
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 2:
        return False, "Missing long carbon chains"
    
    # Count carbons to verify chains are long enough
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:  # Minimum reasonable size for diradylglycerol
        return False, "Carbon chains too short for diradylglycerol"
    
    # Count rotatable bonds to verify flexibility of chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Not enough rotatable bonds for typical chain flexibility"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:  # Minimum reasonable weight for diradylglycerol
        return False, "Molecular weight too low for diradylglycerol"
        
    return True, "Contains glycerol backbone with two substituent groups and one free position"