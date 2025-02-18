"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: CHEBI:xxxxx polypyrrole
"""
from rdkit import Chem
from rdkit.Chem import Mol

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole contains two or more pyrrole units (5-membered conjugated rings with one nitrogen).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Ensure molecule is sanitized to perceive conjugation
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL)
    
    # SMARTS for 5-membered ring with exactly one nitrogen atom
    pyrrole_smarts = '[n]1[c][c][c][c]1'
    pyrrole_pattern = Chem.MolFromSmarts(pyrrole_smarts)
    if not pyrrole_pattern:
        return False, "Failed to initialize pattern"
    
    # Find all matches for 5-membered rings with one nitrogen
    matches = mol.GetSubstructMatches(pyrrole_pattern)
    unique_rings = set()
    conjugated_rings = 0
    
    for match in matches:
        # Check if all bonds in the ring are conjugated
        ring_atoms = set(match)
        conjugated = True
        # Iterate through bonds in the ring
        for i in range(len(match)):
            a1 = match[i]
            a2 = match[(i+1) % len(match)]
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if not bond or not bond.GetIsConjugated():
                conjugated = False
                break
        if conjugated:
            # Ensure unique rings by sorted atom indices
            unique_rings.add(tuple(sorted(match)))
    
    pyrrole_count = len(unique_rings)
    
    if pyrrole_count >= 2:
        return True, f"Contains {pyrrole_count} pyrrole rings"
    else:
        return False, f"Found {pyrrole_count} pyrrole rings, need at least 2"