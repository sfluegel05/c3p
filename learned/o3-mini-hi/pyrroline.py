"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: Organic heteromonocyclic compounds based on a dihydropyrrole (pyrroline) core.
A pyrroline (for our purposes) is defined as a 5‐membered ring that:
  - contains exactly one nitrogen,
  - has exactly one multiple bond that is either:
      • an internal (in–ring) non‐aromatic double bond or
      • an exocyclic double bond from a ring atom to an oxygen or sulfur.
Additionally, to avoid false positives from heavily fused systems we require:
  - if the ring is fused (i.e. not isolated), then the ring should not have an excessively high number of exocyclic substituents.
Because many complex molecules contain fused rings or additional unsaturation, this heuristic is not perfect.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_pyrroline(smiles: str):
    """
    Determines if the given molecule (via its SMILES string) contains a pyrroline (dihydropyrrole)
    core by using the following heuristic:
      1. The molecule is parsed and (when possible) kekulized so that bond orders are explicit.
      2. All 5‐membered rings are retrieved.
      3. For each candidate ring, verify that:
           - It contains exactly one nitrogen.
           - None of its bonds are aromatic (a dihydropyrrole core should be non–aromatic).
           - Count intraring double bonds (explicit double bonds) and also exocyclic double bonds from 
             a ring atom to a heteroatom O or S.
           - The total number of such multiple bonds equals exactly one.
      4. To reduce false positives (especially for fused systems), we compute two metrics:
           a. The number of ring atoms that are unique (appear in only this ring).
           b. The total number of substituents (neighbors outside the ring) attached to ring atoms.
         When the ring is fused (i.e. not all atoms are unique) and has many exocyclic attachments,
         then we reject it.
    If such a ring is found the function returns True plus a message indicating whether the ring is isolated
    or (if fused) acceptable based on its substituent count; otherwise False is returned.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a qualifying pyrroline core is found, False otherwise.
        str: A reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Try to kekulize so explicit bond orders are available.
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass  # If kekulization fails, continue using the molecule as-is.
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found in the molecule"
    
    # Precompute for each atom: in how many rings does it participate.
    atom_ring_counts = {atom.GetIdx(): 0 for atom in mol.GetAtoms()}
    for ring in atom_rings:
        for idx in ring:
            atom_ring_counts[idx] += 1
    
    # Iterate over all 5-membered rings.
    for ring in atom_rings:
        if len(ring) != 5:
            continue
        
        # Count nitrogen atoms in this ring.
        n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_nitrogen != 1:
            continue
        
        # Identify bonds that connect only ring atoms.
        ring_atom_set = set(ring)
        intra_ring_bonds = []
        ring_bonds_aromatic = False
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_atom_set and a2 in ring_atom_set:
                intra_ring_bonds.append(bond)
                if bond.GetIsAromatic():
                    ring_bonds_aromatic = True
        # Skip ring if any bond is aromatic.
        if ring_bonds_aromatic:
            continue
        # If number of bonds found is not 5, something is off.
        if len(intra_ring_bonds) != 5:
            continue
        
        # Count intraring double (non-aromatic) bonds.
        intra_double = sum(1 for bond in intra_ring_bonds 
                           if bond.GetBondType() == rdchem.BondType.DOUBLE)
        
        # Count exocyclic double bonds.
        exo_double = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                # Only consider bond if the other atom is outside the ring.
                other = bond.GetOtherAtom(atom)
                if other.GetIdx() in ring_atom_set:
                    continue
                if bond.GetBondType() == rdchem.BondType.DOUBLE and other.GetAtomicNum() in (8, 16):
                    exo_double += 1
        
        total_multiple = intra_double + exo_double
        if total_multiple != 1:
            continue
        
        # Determine if the ring is isolated or fused.
        unique_count = sum(1 for idx in ring if atom_ring_counts[idx] == 1)
        if unique_count == 5:
            ring_type = "isolated"
        else:
            ring_type = "fused"
        
        # Count exocyclic substituents (neighbors of ring atoms that lie outside the ring).
        exo_substituents = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            exo_substituents += sum(1 for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring_atom_set)
        
        # For fused rings, if there are too many exocyclic attachments,
        # then this ring is likely part of a larger polycyclic system and not a stand–alone dihydropyrrole.
        if ring_type == "fused" and exo_substituents > 4:
            continue
        
        # If we reached here, we found a ring meeting our criteria.
        return True, f"Found a {ring_type} 5‐membered ring with one nitrogen and one multiple bond (dihydropyrrole core)."
    
    return False, "No 5‐membered pyrroline (dihydropyrrole) core meeting our criteria found."

# Example usage:
if __name__ == '__main__':
    # Test using a known pyrroline: 2-acetyl-1-pyrroline
    test_smiles = "CC(=O)C1=NCCC1"
    result, reason = is_pyrroline(test_smiles)
    print(result, reason)