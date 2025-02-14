"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: germacranolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    Germacranolides are sesquiterpene lactones with a 10-membered ring and a lactone.

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

    # 1. Check for the 10-membered ring
    # Look for a 10-membered ring.
    ring10_pattern = Chem.MolFromSmarts("[R10]")
    ring_matches = mol.GetSubstructMatches(ring10_pattern)
    
    if not ring_matches:
        return False, "No 10-membered ring found"
    
    # Ensure only one 10-membered ring exists
    num_10_member_rings = 0
    for match in mol.GetSubstructMatches(ring10_pattern):
      if len(match) == 10:
        num_10_member_rings += 1
    if num_10_member_rings != 1:
        return False, f"Found {num_10_member_rings} 10-membered rings; need exactly 1."

    # 2. Check for the presence of a lactone ring (cyclic ester)
    # The lactone ring could be fused, so we use a generalized SMARTS pattern
    lactone_pattern = Chem.MolFromSmarts("C(=O)O[C;R]") # Check for C=O-O within a ring
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone ring found"
        
    # 3. Molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 350:  # Sesquiterpenes are typically ~C15
        return False, f"Molecular weight ({mol_wt:.2f}) outside typical sesquiterpene range"

    # 4. Ring count
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 2:
        return False, "Must have at least two rings"

    # If all checks pass, then classify as germacranolide
    return True, "Contains a 10-membered ring and a lactone, consistent with a germacranolide"