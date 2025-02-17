"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
#!/usr/bin/env python3
"""
Classifies: 11β-hydroxy steroid
Defined as: Any 11-hydroxy steroid in which the hydroxy group at position 11 has beta- configuration.

This implementation uses a SMARTS pattern to look directly for a beta-oriented hydroxyl group,
i.e. a chiral carbon with beta descriptor (‘@@’) that is directly bonded to an –OH group.
In addition, the molecule must display a fused steroid nucleus (at least 4 rings with at least 3 rings of size 6).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11β-hydroxy steroid based on its SMILES string.
    
    Heuristics:
      1. The molecule must have a fused steroid nucleus (≥4 rings and at least 3 6-membered rings).
      2. It must contain at least one beta-oriented hydroxyl group (i.e. a chiral carbon, marked by @@,
         that is directly bonded to an –OH group) and that carbon is part of one of the 6-membered rings.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria of an 11β-hydroxy steroid, False otherwise.
        str: A message explaining the classification decision.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string."
    except Exception as e:
        return False, f"Error parsing SMILES: {e}"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Check for a steroid-like fused ring system.
    if len(rings) < 4:
        return False, f"Found only {len(rings)} rings; a typical steroid has at least 4 fused rings."
    
    # Count rings of size 6 (steroid nucleus typically has 3 or 4 six-membered rings)
    rings_6 = [ring for ring in rings if len(ring) == 6]
    if len(rings_6) < 3:
        return False, "Fewer than 3 rings of size 6; does not have a typical steroid nucleus."
    
    # SMARTS to find a beta-oriented hydroxy group:
    # [C@@] specifies a chiral carbon with beta configuration.
    # [OX2H] specifies an oxygen with two connections (typically -OH).
    beta_oh_smarts = "[C@@]([OX2H])"
    beta_oh_query = Chem.MolFromSmarts(beta_oh_smarts)
    
    if beta_oh_query is None:
        return False, "Error constructing SMARTS query."
    
    matches = mol.GetSubstructMatches(beta_oh_query)
    if not matches:
        return False, "No beta-oriented hydroxyl group (as defined by SMARTS '[C@@]([OX2H])') found in the molecule."
    
    # Check that at least one of the matching chiral carbons is part of a 6-membered ring.
    for match in matches:
        carbon_idx = match[0]  # first atom in the match is the chiral carbon
        # Check if this carbon is in any 6-membered ring
        for ring in rings_6:
            if carbon_idx in ring:
                return True, "Steroid nucleus detected and at least one beta-oriented hydroxyl group on a 6-membered ring found."
    
    return False, "Matched beta-oriented hydroxyl group not found on any 6-membered ring."

# Example usage and tests.
if __name__ == "__main__":
    test_smiles_list = [
        # prednisolone phosphate
        "C[C@@]12C[C@@H]([C@]3([C@]([C@]2([H])CC[C@@]1(C(COP(=O)(O)O)=O)O)(CCC=4[C@]3(C)C=CC(C4)=O)[H])[H])O",
        # 5beta-dihydroaldosterone
        "C[C@]12CCC(=O)C[C@H]1CC[C@H]1[C@@H]3CC[C@H](C(=O)CO)[C@]3(C[C@H](O)[C@H]21)C=O",
        # withalongolide C
        "C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3[C@@H](O)C[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1"
    ]
    
    for sm in test_smiles_list:
        result, reason = is_11beta_hydroxy_steroid(sm)
        print("SMILES:", sm)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 80)