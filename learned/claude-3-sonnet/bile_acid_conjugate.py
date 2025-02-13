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
        
    # Check for steroid core (four connected rings) - more flexible pattern
    steroid_core = Chem.MolFromSmarts("C1C(C)CCC2C1CCC3C2CC(C)C4CCC(C)C34")
    steroid_core2 = Chem.MolFromSmarts("C1CC2CCC3C(C2C1)CCC4C3(C)CCC4") # Alternative pattern
    
    if not (mol.HasSubstructMatch(steroid_core) or mol.HasSubstructMatch(steroid_core2)):
        return False, "No steroid core structure found"

    # Check for conjugation patterns
    
    # 1. Amino acid conjugations (multiple patterns)
    amino_acid_patterns = [
        Chem.MolFromSmarts("[CX3](=O)[NX3][CX4][CX3](=O)[OX2H1,OX1-]"), # General amino acid
        Chem.MolFromSmarts("[CX3](=O)[NX3]CC(=O)[OX2H1,OX1-]"), # Glycine
        Chem.MolFromSmarts("[CX3](=O)[NX3]C(C)C(=O)[OX2H1,OX1-]"), # Alanine
        Chem.MolFromSmarts("[CX3](=O)[NX3]C(CC(=O)[OX2H1,OX1-])C(=O)[OX2H1,OX1-]"), # Aspartic acid
        Chem.MolFromSmarts("[CX3](=O)[NX3]C(Cc1ccccc1)C(=O)[OX2H1,OX1-]"), # Phenylalanine
    ]
    
    # 2. Taurine conjugation patterns
    taurine_patterns = [
        Chem.MolFromSmarts("[CX3](=O)[NX3]CCS(=O)(=O)[OX2H1,OX1-]"),
        Chem.MolFromSmarts("[CX3](=O)[NX3]CCO[SX4](=O)(=O)[OX2H1,OX1-]")
    ]
    
    # 3. Sulfate conjugation
    sulfate_pattern = Chem.MolFromSmarts("[OX2]S(=O)(=O)[OX2H1,OX1-]")
    
    # 4. Glucuronic acid conjugation
    glucuronic_patterns = [
        Chem.MolFromSmarts("OC1OC(C(=O)O)C(O)C(O)C1O"),
        Chem.MolFromSmarts("OC1OC(C(=O)[OH1,O-])C(O)C(O)C1")
    ]
    
    # 5. Sugar conjugation (more general pattern)
    sugar_pattern = Chem.MolFromSmarts("OC1OC(CO)C(O)C(O)C1")
    
    # Check for conjugations
    found_conjugations = []
    
    # Check amino acid patterns
    for pattern in amino_acid_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            found_conjugations.append("amino acid")
            break
            
    # Check taurine patterns
    for pattern in taurine_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            found_conjugations.append("taurine")
            break
            
    # Check other patterns
    if sulfate_pattern is not None and mol.HasSubstructMatch(sulfate_pattern):
        found_conjugations.append("sulfate")
        
    for pattern in glucuronic_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            found_conjugations.append("glucuronic acid")
            break
            
    if sugar_pattern is not None and mol.HasSubstructMatch(sugar_pattern):
        found_conjugations.append("sugar")
            
    if not found_conjugations:
        return False, "No conjugation patterns found"
        
    # Check for hydroxyl groups (common in bile acids)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_matches < 1:
        return False, "No hydroxyl groups found on steroid core"
    
    # Additional checks for molecular properties
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if mol_weight < 300 or mol_weight > 1200:
        return False, "Molecular weight outside typical range for bile acid conjugates"
    
    return True, f"Bile acid conjugated with {', '.join(set(found_conjugations))}"