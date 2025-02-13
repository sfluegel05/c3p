"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: CHEBI:38083 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
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
    
    # Check for 20 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, "Does not have 20 carbons"
    
    # Check for carbonyl or hydroxyl groups (oxidation)
    has_oxo = any(atom.GetAtomicNum() == 8 and atom.GetExplicitDegree() == 1 for atom in mol.GetAtoms())
    has_oh = any(atom.GetAtomicNum() == 8 and atom.GetExplicitDegree() == 2 for atom in mol.GetAtoms())
    if not (has_oxo or has_oh):
        return False, "No carbonyl or hydroxyl groups found (not oxidized)"
    
    # Check for multiple unsaturations (from EPA, AA, DGLA)
    n_unsaturated = rdMolDescriptors.CalcNumUnsaturatedBonds(mol)
    if n_unsaturated < 3:
        return False, "Not enough unsaturations for EPA/AA/DGLA derivative"
    
    # Look for rings (prostaglandins have rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 0:
        return False, "No rings found (prostaglandin-like structure expected)"
    
    # If all checks pass, classify as icosanoid
    return True, "Contains 20 carbons, oxidized, polyunsaturated, and has rings"