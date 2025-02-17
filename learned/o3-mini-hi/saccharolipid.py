"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: Saccharolipid – Lipids that contain a carbohydrate moiety.
This function uses heuristic SMARTS patterns to detect a sugar ring (as a carbohydrate moiety)
and a long hydrophobic (aliphatic) chain (as part of the lipid backbone).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is defined as a lipid that contains a carbohydrate moiety.
    Heuristics used:
      - A carbohydrate moiety is approximated by the presence of a pyranose or furanose ring,
        i.e. a 6-membered ring (pyranose) or 5-membered ring (furanose) containing oxygen(s)
        and several hydroxyl groups.
      - A lipid is approximated by the presence of a long contiguous aliphatic chain (8 or more carbon atoms).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a saccharolipid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for a pyranose-like ring (6-membered sugar ring).
    # This pattern is a heuristic—it seeks a ring with 5 carbons and 1 oxygen in a typical sugar ring,
    # with several hydroxyl substituents.
    pyranose_smarts = "[OX2r6][C@H]1[C@@H]([OX2H])[C@H]([OX2H])[C@H]([OX2H])[C@H]1[OX2H]"
    pyranose_pattern = Chem.MolFromSmarts(pyranose_smarts)
    
    # SMARTS pattern for a furanose-like ring (5-membered sugar ring).
    furanose_smarts = "[OX2r5][C@H]1[C@@H]([OX2H])[C@H]([OX2H])[C@H]1[OX2H]"
    furanose_pattern = Chem.MolFromSmarts(furanose_smarts)
    
    sugar_found = False
    # Try to find a pyranose-like sugar ring.
    if mol.HasSubstructMatch(pyranose_pattern):
        sugar_found = True
    # If not found, try a furanose-like sugar ring.
    elif mol.HasSubstructMatch(furanose_pattern):
        sugar_found = True

    if not sugar_found:
        return False, "No carbohydrate (sugar ring) moiety found"
    
    # Define a simple SMARTS to detect a long aliphatic chain.
    # Here we check for a chain of at least 8 contiguous carbon atoms.
    long_chain_smarts = "CCCCCCCC"
    long_chain_pattern = Chem.MolFromSmarts(long_chain_smarts)
    
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long aliphatic (lipid-like) chain found"
    
    # Optionally, one could perform further checks (e.g. molecular weight or rotatable bonds)
    # to ensure the lipid part is sizeable.
    
    return True, "Molecule has both a carbohydrate moiety and a lipid-like aliphatic chain, classifying it as a saccharolipid"

# Below are some example calls (which you might uncomment to test if needed).
# For example, one could test one of the provided saccharolipid SMILES.
# example_smiles = "[H][C@@]1(O[C@@](C[C@@H](O)[C@H]1O)..."
# print(is_saccharolipid(example_smiles))