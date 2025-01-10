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
        
    # Check for uracil ring (more flexible pattern)
    uracil_pattern = Chem.MolFromSmarts("[#7]1[#6]=[#6][#6](=[O])[#7][#6]1=[O]")
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "No uracil ring found"
    
    # Check for ribose connected to uracil (less strict stereochemistry)
    ribose_pattern = Chem.MolFromSmarts("[OH1,O][CH2]1O[CH]([CH]([OH1,O])[CH]1[OH1,O])n1ccc(=O)[nH]c1=O")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose moiety found connected to uracil"
        
    # Check for diphosphate linkage (more flexible pattern)
    diphosphate_pattern = Chem.MolFromSmarts("[O,OH]-[P](=O)([O,OH])-O-[P](=O)([O,OH])-[O,OH]")
    if not mol.HasSubstructMatches(diphosphate_pattern):
        return False, "No diphosphate linkage found"
    
    # Check for complete UDP core with proper connectivity
    udp_core_pattern = Chem.MolFromSmarts("""
        [CH2]1O[CH]([CH]([OH1,O])[CH]1[OH1,O])n1ccc(=O)[nH]c1=O  # Ribose-uracil
        .[CH2]OP(O)(=O)OP(O)(=O)O                                 # Phosphate connection
    """)
    if not mol.HasSubstructMatch(udp_core_pattern):
        return False, "Missing or incorrect UDP core structure"
    
    # Check for sugar moiety connected via diphosphate
    # More specific sugar pattern that matches common sugar structures
    sugar_connection_pattern = Chem.MolFromSmarts("""
        [CH2]1O[CH]([CH]([OH1,O])[CH]1[OH1,O])n1ccc(=O)[nH]c1=O  # UDP part
        .[CH2]OP(O)(=O)OP(O)(=O)O[CH]1[O,N][CH]([CH,NH,NC])([CH,OH,O])[CH,C]([CH,OH,O,N])[CH,C]1  # Sugar connection
    """)
    if not mol.HasSubstructMatch(sugar_connection_pattern):
        return False, "No sugar moiety found connected via diphosphate"

    # Additional check for common sugar modifications
    sugar_modifications = [
        (Chem.MolFromSmarts("C(=O)[OH]"), "carboxylic acid"),  # Uronic acids
        (Chem.MolFromSmarts("NC(=O)C"), "N-acetyl"),           # N-acetyl sugars
        (Chem.MolFromSmarts("S(=O)(=O)[OH]"), "sulfate"),      # Sulfated sugars
        (Chem.MolFromSmarts("[NH2]"), "amino"),                # Amino sugars
    ]
    
    modifications = []
    for pattern, mod_type in sugar_modifications:
        if pattern and mol.HasSubstructMatch(pattern):
            modifications.append(mod_type)
            
    base_message = "Contains UDP moiety connected to sugar via diphosphate linkage"
    if modifications:
        return True, f"{base_message} with {', '.join(modifications)} modification(s)"
    return True, base_message