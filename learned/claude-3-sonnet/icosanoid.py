"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: CHEBI:38083 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    An icosanoid is a signaling molecule derived from oxidation of C20 essential fatty acids
    like icosapentaenoic acid (EPA), arachidonic acid (AA), and dihomo-gamma-linolenic acid (DGLA).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for 20 or fewer carbons (some examples have less than 20)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count > 20:
        return False, "Has more than 20 carbons"
    
    # Check for carbonyl or hydroxyl groups (oxidation)
    oxo_pattern = Chem.MolFromSmarts("[C=O]")
    has_oxo = mol.HasSubstructMatch(oxo_pattern)
    oh_pattern = Chem.MolFromSmarts("[OX1H]")
    has_oh = mol.HasSubstructMatch(oh_pattern)
    if not (has_oxo or has_oh):
        return False, "No carbonyl or hydroxyl groups found (not oxidized)"
    
    # Check for multiple unsaturations (from EPA, AA, DGLA)
    n_unsaturated = rdMolDescriptors.CalcNumUnsaturatedBonds(mol)
    if n_unsaturated < 2:
        return False, "Not enough unsaturations for EPA/AA/DGLA derivative"
    
    # Look for rings (many, but not all, examples have rings)
    ring_info = mol.GetRingInfo()
    has_rings = ring_info.NumRings() > 0
    
    # If all checks pass, classify as icosanoid
    if has_rings:
        reason = "Contains oxidized, polyunsaturated C20 (or fewer) chain with rings"
    else:
        reason = "Contains oxidized, polyunsaturated C20 (or fewer) chain"
    return True, reason