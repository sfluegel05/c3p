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
    
    # Basic flavan-3-ol core structure (2-phenyl-3,4-dihydro-2H-chromene-3-ol)
    # Note: The SMARTS pattern matches the core structure with potential substitutions
    flavan_3_ol_pattern = Chem.MolFromSmarts(
        "[OX2H0,OX2H1]-1-[#6]-2-[#6]-[#6]-[#6](=[#6]-[#6](=[#6]-2)-[OX2H0,OX2H1])-[#6]-1-[#6]-1=[#6]-[#6]=[#6]-[#6]=[#6]-1"
    )
    
    if not mol.HasSubstructMatch(flavan_3_ol_pattern):
        return False, "No flavan-3-ol core structure found"
    
    # Check for hydroxyl group at C3 position (characteristic of flavan-3-ols)
    c3_oh_pattern = Chem.MolFromSmarts("[OX2H]-[CH1X4]-1-[#6]-2-[#6]~[#6]~[#6]~[#6]~[#6]-2-O-[CH1X4]-1-[#6]")
    if not mol.HasSubstructMatch(c3_oh_pattern):
        return False, "Missing hydroxyl group at C3 position"
    
    # Count hydroxyl groups - catechins typically have multiple OH groups
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 3:
        return False, f"Too few hydroxyl groups ({oh_matches}) for a catechin"
    
    # Look for common substitution patterns
    
    # Check for gallate ester (common in many catechins)
    gallate_pattern = Chem.MolFromSmarts("[OX2]-C(=O)-c1c(O)c(O)c(O)cc1")
    has_gallate = mol.HasSubstructMatch(gallate_pattern)
    
    # Check for typical A-ring hydroxylation pattern (5,7-dihydroxy)
    a_ring_pattern = Chem.MolFromSmarts("O-1-c2c(O)cc(O)cc2CC[CH1]-1")
    has_typical_a_ring = mol.HasSubstructMatch(a_ring_pattern)
    
    # Check for B-ring catechol pattern (3',4'-dihydroxy)
    b_ring_pattern = Chem.MolFromSmarts("c1c(O)c(O)ccc1")
    has_typical_b_ring = mol.HasSubstructMatch(b_ring_pattern)
    
    # Additional characteristics
    features = []
    if has_gallate:
        features.append("gallate ester")
    if has_typical_a_ring:
        features.append("5,7-dihydroxy pattern")
    if has_typical_b_ring:
        features.append("3',4'-dihydroxy pattern")
    
    # Verify molecular weight is in reasonable range for catechins (290-1000 Da)
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 1000:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for catechins"
    
    feature_str = ", ".join(features) if features else "basic"
    return True, f"Contains flavan-3-ol core with {oh_matches} hydroxyl groups. Type: {feature_str} catechin"