"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_bile_acid_conjugate, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Multiple patterns for steroid core to catch different variations
    steroid_patterns = [
        # Basic steroid core with flexible ring connections
        "[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]3~[#6]~[#6]2~[#6]~1",
        # Alternative pattern allowing for keto groups
        "[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6,#8]~[#6]~[#6]~[#6]~4~[#6]~[#6]3~[#6]~[#6]2~[#6]~1",
        # More specific bile acid core pattern
        "[#6]1~[#6]~[#6]~[#6]2~[#6]([#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]3~[#6]~2)~[#6]~1"
    ]
    
    found_core = False
    for pattern in steroid_patterns:
        core_pattern = Chem.MolFromSmarts(pattern)
        if core_pattern is not None and mol.HasSubstructMatch(core_pattern):
            found_core = True
            break
            
    if not found_core:
        return False, "No steroid core structure found"

    # Conjugation patterns
    conjugation_patterns = {
        "glycine": ["[CX3](=O)[NX3]CC(=O)[OX2H1,OX1-]"],
        "taurine": ["[CX3](=O)[NX3]CCS(=O)(=O)[OX2H1,OX1-]",
                   "[CX3](=O)[NX3]CCO[SX4](=O)(=O)[OX2H1,OX1-]"],
        "amino acid": ["[CX3](=O)[NX3][CX4][CX3](=O)[OX2H1,OX1-]",
                      "[CX3](=O)[NX3]C([CX4])[CX3](=O)[OX2H1,OX1-]"],
        "sulfate": ["[OX2]S(=O)(=O)[OX2H1,OX1-]"],
        "glucuronic acid": ["OC1OC(C(=O)O)C(O)C(O)C1",
                           "OC1OC(C(=O)[OH1,O-])C(O)C(O)C1O"],
        "sugar": ["OC1OC(CO)C(O)C(O)C1"]
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
    
    # Molecular weight check
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if mol_weight < 300 or mol_weight > 1200:
        return False, "Molecular weight outside typical range for bile acid conjugates"
    
    # Additional check for carboxylic acid group (common in bile acids)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")
    if carboxyl_pattern is not None and not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    return True, f"Bile acid conjugated with {', '.join(set(found_conjugations))}"