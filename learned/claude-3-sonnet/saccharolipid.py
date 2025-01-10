"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: saccharolipids - lipids containing carbohydrate moieties
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    Saccharolipids are lipids that contain carbohydrate moieties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for minimum complexity - saccharolipids are large molecules
    if mol.GetNumAtoms() < 20:
        return False, "Molecule too small to be a saccharolipid"

    # Look for sugar rings (pyranose/furanose)
    sugar_pattern = Chem.MolFromSmarts("[CR1]1[CR1][CR1][CR1][CR1][OR1]1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar rings found"

    # Look for long carbon chains (lipid part)
    lipid_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2]")
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if not lipid_matches:
        return False, "No long carbon chains found"

    # Count key elements typical for saccharolipids
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if num_carbons < 12:
        return False, "Insufficient carbon atoms for saccharolipid"
    if num_oxygens < 6:
        return False, "Insufficient oxygen atoms for saccharolipid"

    # Look for ester/amide linkages (common in saccharolipids)
    ester_pattern = Chem.MolFromSmarts("[#6]-C(=O)-O-[#6]")
    amide_pattern = Chem.MolFromSmarts("[#6]-C(=O)-N-[#6]")
    
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_amide = mol.HasSubstructMatch(amide_pattern)
    
    if not (has_ester or has_amide):
        return False, "No ester or amide linkages found"

    # Check for hydroxyl groups (characteristic of sugars)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Insufficient hydroxyl groups for sugar moiety"

    # Optional: Check for phosphate groups (common in many saccharolipids)
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)([O-,OH])([O-,OH])")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)

    # Calculate molecular weight - saccharolipids are typically large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for saccharolipid"

    # Construct reason string
    reason = f"Contains {len(sugar_matches)} sugar ring(s), lipid chains, "
    reason += f"{len(hydroxyl_matches)} hydroxyl groups, "
    reason += "and appropriate linkages between sugar and lipid components"
    if has_phosphate:
        reason += " with phosphate groups"

    return True, reason