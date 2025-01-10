"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: sterol ester
A steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid 
with the 3-hydroxy group of a sterol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core pattern - four fused rings (three 6-membered, one 5-membered)
    steroid_core = Chem.MolFromSmarts("[cR2,CR2]1[cR2,CR2][cR2,CR2][cR2,CR2]2[cR2,CR2][cR2,CR2][cR2,CR2]3[cR2,CR2][cR2,CR2][cR2,CR2]4[cR2,CR2][cR2,CR2][cR2,CR2][cR2,CR2]4[cR2,CR2]3[cR2,CR2]2[cR2,CR2]1")
    if steroid_core is None:
        return False, "Invalid steroid core SMARTS pattern"
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    if ester_pattern is None:
        return False, "Invalid ester SMARTS pattern"
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"

    # Pattern for 3-position ester on steroid core
    sterol_3_ester = Chem.MolFromSmarts("[C,c]1[C,c][C,c][C,c]([OX2][CX3](=[OX1])[#6])[C,c]2")
    if sterol_3_ester is None:
        return False, "Invalid sterol ester SMARTS pattern"
    if not mol.HasSubstructMatch(sterol_3_ester):
        return False, "Ester not at position 3 of steroid core"

    # Look for characteristic sterol side chain (typically 8 or more carbons)
    side_chain = Chem.MolFromSmarts("[CH2][CH2][CH][CH]([CH3])[CH2][CH2][CH3]")
    if side_chain is None:
        return False, "Invalid side chain SMARTS pattern"
    if not mol.HasSubstructMatch(side_chain):
        return False, "Missing characteristic sterol side chain"

    # Count carbons (sterols typically have 27-29 carbons plus ester group)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27:
        return False, f"Too few carbons ({c_count}) for a sterol ester"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Too few rings for a sterol structure"

    # Look for angular methyl groups characteristic of sterols
    angular_methyl = Chem.MolFromSmarts("[CH3][C]([C,c])([C,c])[C,c]")
    if angular_methyl is None:
        return False, "Invalid angular methyl SMARTS pattern"
    methyl_matches = mol.GetSubstructMatches(angular_methyl)
    if len(methyl_matches) < 2:
        return False, "Missing characteristic angular methyl groups"

    return True, "Contains sterol core with ester at position 3 and characteristic side chain"