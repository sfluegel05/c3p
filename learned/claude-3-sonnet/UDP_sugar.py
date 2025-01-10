"""
Classifies: CHEBI:17297 UDP-sugar
"""
"""
Classifies: UDP-sugar
A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached 
to an unspecified sugar via an anomeric diphosphate linkage.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for uracil ring with specific pattern found in UDP sugars
    uracil_pattern = Chem.MolFromSmarts("n1ccc(=O)[nH]c1=O")
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "No uracil ring found"
    
    # Check for ribose attached to uracil
    # Pattern: ribose with hydroxyl groups and connection to uracil
    uridine_pattern = Chem.MolFromSmarts("[OH1][C@H]1O[C@H]([C@H]([OH1])[C@@H]1[OH1])n1ccc(=O)[nH]c1=O")
    if not mol.HasSubstructMatch(uridine_pattern):
        return False, "No uridine moiety found"
        
    # Check for diphosphate group with specific UDP pattern
    diphosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)O")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No diphosphate linkage found"
        
    # Check for sugar moiety - more flexible to account for various sugars
    sugar_pattern = Chem.MolFromSmarts("[C][O,N][C]([!$(C=O)])[C]([OH,NH2,NC,O])[C]")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moiety found"
        
    # Verify complete UDP core structure with correct connectivity
    udp_core = Chem.MolFromSmarts("""
        [C][C]1O[C@H]([CH2]OP(O)(=O)OP(O)(=O)O)[C@H]([OH1])[C@@H]1[OH1]  # Ribose-phosphate
        .[#7]1[#6](=[O])[#7][#6][#6]1=[O]                                 # Uracil
    """)
    if not mol.HasSubstructMatch(udp_core):
        return False, "Missing or incorrect UDP core structure"
    
    # Check for proper connection between UDP and sugar
    connection_pattern = Chem.MolFromSmarts("[C]OP(O)(=O)OP(O)(=O)O[C]")
    if not mol.HasSubstructMatch(connection_pattern):
        return False, "UDP not properly connected to sugar"
    
    # Additional check for anomeric linkage to sugar
    # The sugar should be connected via an anomeric carbon
    anomeric_pattern = Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)O[C]1[O][C]")
    if not mol.HasSubstructMatch(anomeric_pattern):
        return False, "Sugar not connected via anomeric position"

    return True, "Contains UDP moiety connected to sugar via anomeric diphosphate linkage"