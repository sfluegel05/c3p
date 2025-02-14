"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone with a germacrane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define broader SMARTS pattern for germacrane skeleton
    germacrane_patterns = [
        Chem.MolFromSmarts("C1CCC(CC1)C2CCC(CC2)C3=CC=CC=C3"), # Flexible pattern
        Chem.MolFromSmarts("C1CCC2CCCC2C1")  # Decalin or other common motifs
    ]
    
    # Attempt to match any germacrane pattern
    germacrane_match = any(mol.HasSubstructMatch(pattern) for pattern in germacrane_patterns)
    if not germacrane_match:
        return False, "No germacrane-like structure found"

    # Ensure it contains a lactone group in a 10-membered ring (common feature)
    lactone_pattern = Chem.MolFromSmarts("C1OC(=O)C(C2CCCCC2)=C1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No appropriate or ring-contained lactone group found"
    
    # General ring checks for sesquiterpene-like properties
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized() or ring_info.NumRings() < 2:
        return False, "Too few rings for complex germacrane structure"

    # Consider oxygen-containing functional groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "Too few oxygens for typical germacranolide functionalization"

    # General complexity or molecular weight check (optional, depending on data constraints)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for typical germacranolide"
    
    return True, "Contains germacrane skeleton and key functional groups typical of germacranolides"