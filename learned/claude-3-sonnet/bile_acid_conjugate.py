"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate consists of a bile acid core (steroid structure) with
    conjugating groups like amino acids, taurine, sulfate, or sugars.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_bile_acid_conjugate, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Bile acid core patterns - more specific to common bile acid structures
    steroid_patterns = [
        # Basic steroid core (more flexible)
        "[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]3~[#6]~[#6]2~[#6]~1",
        # 5β-cholane core specific pattern
        "[CH2]1[CH2][CH2][C@H]2[CH2][CH2][C@H]3[CH2][CH2][C@]4([CH2][CH2][CH2][CH2]4)[C@H]3[CH2][C@H]2[CH2]1",
        # More general pattern allowing for oxidized positions
        "[#6]1[#6][#6][#6]2[#6][#6][#6]3[#6][#6][#6]4[#6,#8][#6][#6][#6]4[#6][#6]3[#6][#6]2[#6]1",
        # Pattern allowing for common substituents
        "[#6]1[#6][#6][#6]2[#6]([#6][#6]3[#6][#6][#6]4[#6][#6][#6,#8][#6]4[#6]3[#6]2)[#6]1",
        # Specific pattern for common bile acid oxidation states
        "[CH2]1[CH2][CH2][CH]2[CH2][CH2][C]3([CH2][CH2][C]4([CH2][CH2][CH2][CH2]4)[C]3[CH2][C]2[CH2]1)C"
    ]
    
    found_core = False
    for pattern in steroid_patterns:
        core_pattern = Chem.MolFromSmarts(pattern)
        if core_pattern is not None and mol.HasSubstructMatch(core_pattern):
            found_core = True
            break
            
    if not found_core:
        return False, "No bile acid core structure found"

    # Conjugation patterns
    conjugation_patterns = {
        "glycine": ["NCC(=O)[OH]", "[CX3](=O)[NX3]CC(=O)[OX2H1,OX1-]"],
        "taurine": ["[CX3](=O)[NX3]CCS(=O)(=O)[OX2H1,OX1-]",
                   "NCCS(=O)(=O)[OH]"],
        "amino acid": ["[CX3](=O)[NX3][CH]([CX4])C(=O)[OX2H1,OX1-]",
                      "[CX3](=O)[NX3]C([CX4])[CX3](=O)[OX2H1,OX1-]",
                      "NC([CX4])C(=O)[OH]"],
        "sulfate": ["OS(=O)(=O)[OH]", "[OX2]S(=O)(=O)[OX2H1,OX1-]"],
        "glucuronic acid": ["[OH]C1O[CH]([CH]([OH])[CH]([OH])[CH]([OH])C1=O)",
                           "OC1OC(C(=O)O)C(O)C(O)C1"],
        "sugar": ["[OH]C1O[CH]([CH]([OH])[CH]([OH])[CH]([OH])CO)"]
    }
    
    found_conjugations = []
    for conj_type, patterns in conjugation_patterns.items():
        for pattern in patterns:
            patt = Chem.MolFromSmarts(pattern)
            if patt is not None and mol.HasSubstructMatch(patt):
                found_conjugations.append(conj_type)
                break
                
    if not found_conjugations:
        return False, "No conjugation patterns found"

    # Check for hydroxyl groups (common in bile acids)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H1]")
    if hydroxyl_pattern is not None:
        hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
        if hydroxyl_matches < 1:
            return False, "No hydroxyl groups found on steroid core"
    
    # Molecular weight check - typical range for bile acid conjugates
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if mol_weight < 300 or mol_weight > 1200:
        return False, "Molecular weight outside typical range for bile acid conjugates"
    
    # Check for carboxylic acid group or conjugated derivative
    carboxyl_pattern = Chem.MolFromSmarts("[$([CX3](=O)[OX2H1]),$([CX3](=O)[OX1-]),$([CX3](=O)[NX3])]")
    if carboxyl_pattern is not None and not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid or conjugated derivative found"
    
    return True, f"Bile acid conjugated with {', '.join(set(found_conjugations))}"