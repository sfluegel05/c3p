"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are flavan-3-ols with various substitution patterns.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic chroman (benzopyran) core - more flexible pattern
    chroman_pattern = Chem.MolFromSmarts("O1CCc2ccccc2C1")
    if not mol.HasSubstructMatch(chroman_pattern):
        return False, "No benzopyran core structure found"

    # Check for the basic flavan core (more flexible pattern)
    # This matches both catechin and epicatechin configurations
    flavan_pattern = Chem.MolFromSmarts("O1[CH1][CH1]c2ccccc2[CH1]1c1ccccc1")
    if not mol.HasSubstructMatch(flavan_pattern):
        return False, "No flavan core structure found"

    # Check for hydroxyl at C3 position (more flexible pattern)
    # Matches both axial and equatorial OH
    c3_oh_pattern = Chem.MolFromSmarts("O1[CH1][CH1](O)c2ccccc2[CH1]1")
    if not mol.HasSubstructMatch(c3_oh_pattern):
        return False, "Missing hydroxyl group at C3 position"

    # Count hydroxyl groups
    oh_pattern = Chem.MolFromSmarts("[OH1]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 2:  # Most catechins have at least 2 OH groups
        return False, f"Too few hydroxyl groups ({oh_matches}) for a catechin"

    # Look for common substitution patterns, but don't require them
    features = []

    # Check for common A-ring substitution (5,7-dihydroxy pattern)
    a_ring_pattern = Chem.MolFromSmarts("O1[CH1][CH1]c2c(O)cc(O)cc2[CH1]1")
    if mol.HasSubstructMatch(a_ring_pattern):
        features.append("5,7-dihydroxy pattern")

    # Check for B-ring patterns (various possible hydroxylation patterns)
    catechol_pattern = Chem.MolFromSmarts("c1c(O)c(O)ccc1")
    pyrogallol_pattern = Chem.MolFromSmarts("c1c(O)c(O)c(O)cc1")
    
    if mol.HasSubstructMatch(pyrogallol_pattern):
        features.append("3',4',5'-trihydroxy (gallocatechin-type)")
    elif mol.HasSubstructMatch(catechol_pattern):
        features.append("3',4'-dihydroxy (catechin-type)")

    # Check for ester derivatives (gallate, coumarate, etc.)
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    if mol.HasSubstructMatch(ester_pattern):
        features.append("ester derivative")

    # Molecular weight check - more permissive range
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1200:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for catechins"

    # Ring count check - catechins should have at least 2 rings
    ring_count = Chem.rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2:
        return False, "Insufficient ring count for a catechin structure"

    feature_str = ", ".join(features) if features else "basic"
    return True, f"Contains flavan-3-ol core with {oh_matches} hydroxyl groups. Features: {feature_str}"