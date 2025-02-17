"""
Classifies: CHEBI:33913 corrinoid
"""
"""
Classifies: a corrinoid
Definition: A derivative of the corrin nucleus, which contains four reduced or partly reduced pyrrole rings
joined in a macrocycle by three =C- groups and one direct carbon-carbon bond linking alpha positions.
This implementation uses heuristic checks.
"""

from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    
    The heuristic used is as follows:
      1. Look for a “macrocycle” (a ring system found by RDKit’s SSSR)
         whose size is between 15 and 25 atoms and that contains exactly 4 nitrogen atoms.
      2. If that fails, look for at least four 5-membered rings (“pyrrole‐like rings”)
         each having exactly one nitrogen atom.
      3. Also if the molecule contains cobalt atoms (common in many corrinoids), this is used
         as supporting evidence.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is judged as a corrinoid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # First heuristic: search for a macrocycle of size 15-25 that contains exactly 4 nitrogens.
    for ring in rings:
        ring_size = len(ring)
        if 15 <= ring_size <= 25:
            n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogen == 4:
                return True, f"Found a macrocycle of size {ring_size} with exactly 4 ring nitrogens."
    
    # Second heuristic: count rings that look like a pyrrole (5-membered ring with one nitrogen).
    pyrrole_like_rings = []
    for ring in rings:
        if len(ring) == 5:
            n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogen == 1:
                pyrrole_like_rings.append(ring)
                
    if len(pyrrole_like_rings) >= 4:
        return True, f"Found {len(pyrrole_like_rings)} pyrrole-like rings (5-membered rings with one nitrogen each), consistent with a corrinoid."
    
    # Bonus check: many corrinoid structures bind a cobalt ion.
    cobalt_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 27]
    if cobalt_atoms and len(pyrrole_like_rings) >= 3:
        return True, "Cobalt present and several pyrrole-like rings detected, consistent with many corrinoids."
    
    return False, "Could not detect a corrin-like macrocycle or sufficient pyrrole-like rings."
    
# Example usage (uncomment to test):
# smiles_example = "C1C\\C2=C\\C3=N\\C(CC3)=C/C3=N/C(CC3)=C\\C3=NC(CC3)C1N2"  # a corrin example
# flag, reason = is_corrinoid(smiles_example)
# print(flag, reason)