"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
#!/usr/bin/env python3
"""
Classifies: 11β-hydroxy steroid
Defined as: Any 11-hydroxy steroid in which the hydroxy group at position 11 has beta- configuration.
Algorithm:
  1. Parse the SMILES string.
  2. Verify that the molecule has a steroid-like fused ring system (at least 4 rings, with at least 3 rings of size 6).
  3. Look for an -OH group (oxygen with one heavy neighbor) where the attached carbon is chiral and
     is a member of a 6-membered ring. Instead of checking if the fragment SMILES starts with "[C@@",
     we now check if it contains the sequence "@@" which is a heuristic for beta orientation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11β-hydroxy steroid based on its SMILES string.
    The heuristic is:
      - The molecule must display a fused steroid nucleus (≥4 rings and at least 3 rings of size 6).
      - There must be at least one hydroxyl (-OH) group attached to a chiral carbon that belongs
        to at least one 6-membered ring. The local stereochemical SMILES for that carbon must include
        the marker "@@" indicating beta orientation.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule passes the heuristic for an 11β-hydroxy steroid.
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
    
    # Count rings of size 6 (as most steroid core rings are 6-membered)
    rings_6 = [ring for ring in rings if len(ring) == 6]
    if len(rings_6) < 3:
        return False, "Fewer than 3 rings of size 6; does not appear to have a typical steroid nucleus."
    
    # Look for a candidate beta-oriented hydroxyl.
    beta_hydroxy_found = False
    # Iterate over oxygen atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        neighbors = atom.GetNeighbors()
        # Check for an -OH group: oxygen bonded to exactly one heavy (non-hydrogen) atom.
        # (We assume implicit H's; if a hydrogen were explicit it would be added to the molecule.)
        if len(neighbors) != 1:
            continue
        heavy = neighbors[0]
        # Verify the bond is a single bond.
        bond = mol.GetBondBetweenAtoms(heavy.GetIdx(), atom.GetIdx())
        if bond is None or bond.GetBondTypeAsDouble() != 1:
            continue
        
        # Check that the attached heavy atom is a carbon.
        if heavy.GetAtomicNum() != 6:
            continue
        
        # Check that this carbon is a chiral center.
        if heavy.GetChiralTag() not in (Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
                                        Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW):
            continue

        # Verify the carbon is part of at least one 6-membered ring
        in_6_ring = False
        for ring in rings:
            if len(ring) == 6 and heavy.GetIdx() in ring:
                in_6_ring = True
                break
        if not in_6_ring:
            continue
        
        # Generate a SMILES fragment for the carbon (with chirality)
        try:
            frag = Chem.MolFragmentToSmiles(mol, atomsToUse=[heavy.GetIdx()], canonical=True, isomericSmiles=True)
        except Exception:
            continue
        
        # Instead of checking that the fragment starts with "[C@@", we now check if it contains "@@"
        # (a heuristic for beta orientation in many steroid SMILES conventions).
        if "@@" in frag:
            beta_hydroxy_found = True
            break

    if not beta_hydroxy_found:
        return False, "No beta-oriented hydroxyl group on a chiral carbon in a 6-membered steroid ring was found."
    
    return True, "Steroid nucleus detected with at least 4 fused rings (≥3 of size 6) and a beta-oriented hydroxyl group likely at C-11."

# Example usage and tests:
if __name__ == "__main__":
    test_smiles_list = [
        # prednisolone phosphate:
        "C[C@@]12C[C@@H]([C@]3([C@]([C@]2([H])CC[C@@]1(C(COP(=O)(O)O)=O)O)(CCC=4[C@]3(C)C=CC(C4)=O)[H])[H])O",
        # 5beta-dihydroaldosterone:
        "C[C@]12CCC(=O)C[C@H]1CC[C@H]1[C@@H]3CC[C@H](C(=O)CO)[C@]3(C[C@H](O)[C@H]21)C=O",
        # withalongolide C:
        "C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3[C@@H](O)C[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1",
    ]
    
    for sm in test_smiles_list:
        result, reason = is_11beta_hydroxy_steroid(sm)
        print("SMILES:", sm)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 80)