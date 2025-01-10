"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid typically has a 15-carbon skeleton derived from three isoprene units,
    often arranged in complex ring structures and functionalized with components like alcohols, ketones, or esters.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check carbon count; sesquiterpenoids have ~15 carbons, allow some variance
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (13 <= c_count <= 18):
        return False, f"Carbon count is {c_count}, expected closer to 15 for sesquiterpenoid"
    
    # Look for complex bicyclic or tricyclic structures
    bicyclic_patterns = [
        Chem.MolFromSmarts("C1CCC2CCC1C2"), #Example generic bicyclic
        Chem.MolFromSmarts("C1C2CCC1CC2"),  #Another possible bicyclic skeleton
        Chem.MolFromSmarts("C1CCC2CC3CCC(C3)C2C1"), #Example tricyclic pattern
    ]
    has_bicyclic = any(mol.HasSubstructMatch(pattern) for pattern in bicyclic_patterns)
    if not has_bicyclic:
        return False, "No complex bicyclic or tricyclic-like structure found"
    
    # Presence of functional groups typical in terpenoids
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OH]")
    ketone_pattern = Chem.MolFromSmarts("C=O")
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    
    if not (mol.HasSubstructMatch(alcohol_pattern) 
            or mol.HasSubstructMatch(ketone_pattern) 
            or mol.HasSubstructMatch(ester_pattern)
            or mol.HasSubstructMatch(carboxylic_pattern)):
        return False, "Missing typical functional groups for sesquiterpenoids like alcohols, ketones, esters, or carboxylic acids"
    
    return True, "Likely contains a 15-carbon sesquiterpenoid skeleton with complex ring structures and functional groups"