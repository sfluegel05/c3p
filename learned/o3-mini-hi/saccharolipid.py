"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: Saccharolipid – Lipids that contain a carbohydrate moiety.
This function uses heuristic ring analysis to detect a sugar ring (carbohydrate moiety)
and a substructure search for a long aliphatic chain for the lipid portion.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is defined as a lipid that contains a carbohydrate moiety.
    
    Heuristics used:
      - A carbohydrate moiety is approximated by detecting at least one
        ring of size 5 (furanose-like) or 6 (pyranose-like) that contains at least one oxygen atom.
      - A lipid is approximated by the presence of a long contiguous aliphatic chain,
        here defined as a chain of 8 or more carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a saccharolipid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a sugar ring using ring info.
    # A sugar ring is typically 5 or 6 atoms; we require that at least one heteroatom is oxygen.
    ring_info = mol.GetRingInfo()
    sugar_found = False
    for ring in ring_info.AtomRings():
        if len(ring) in (5, 6):
            # Count oxygen atoms in this ring.
            oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            # For a typical sugar ring, at least one oxygen should be present.
            if oxygen_count >= 1:
                sugar_found = True
                break
    if not sugar_found:
        return False, "No carbohydrate (sugar ring) moiety found"
    
    # Use a simple SMARTS pattern to detect a long aliphatic chain.
    # For our purposes, a substructure of 8 contiguous carbons is a crude indicator of a lipid chain.
    long_chain_smarts = "CCCCCCCC"
    long_chain_pattern = Chem.MolFromSmarts(long_chain_smarts)
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long aliphatic (lipid-like) chain found"
    
    # If both conditions are met, classify as saccharolipid.
    return True, "Molecule has a sugar ring and a long aliphatic chain, classifying it as a saccharolipid"

# Example usage:
# test_smiles = "CCCCCCCCCCCCCCCC[C@H](C)C[C@H](C)OC1OC(CO)C(O)C(O)C1O"  # A hypothetical saccharolipid SMILES
# print(is_saccharolipid(test_smiles))