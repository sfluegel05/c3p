"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: Organic heteromonocyclic compounds based on a dihydropyrrole (pyrroline) core.
A pyrroline (here) is defined as a 5‐membered ring that has exactly one nitrogen atom and
exactly one multiple bond in the ring. The multiple bond may be an internal (in‐ring) double
bond or an exocyclic double bond (e.g. a C=O or C=S where the C is in the ring).
Note: Because many complex molecules may contain fused rings or exocyclic unsaturation,
this heuristic may over– or under–classify some structures.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_pyrroline(smiles: str):
    """
    Determines if the given molecule (via its SMILES string) contains a pyrroline (dihydropyrrole)
    core. For our purpose the ring is recognized if:
      - It is a 5-membered ring (from the SSSR)
      - It contains exactly one nitrogen atom (in the ring)
      - It has exactly one “multiple bond” associated with the ring.
        This multiple bond is counted as either:
          (a) an intraring bond that is drawn as a double bond
          (b) or an exocyclic double bond from a ring atom to an external oxygen or sulfur.
          
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if a qualifying pyrroline core is found, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Attempt to Kekulize so that double bonds are explicit
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        # Some molecules may fail kekulization. Continue with original structure.
        pass
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # list of tuples of atom indices
    
    if not atom_rings:
        return False, "No rings found in the molecule"
    
    # Loop over each 5-membered ring
    for ring in atom_rings:
        if len(ring) != 5:
            continue  # only consider 5-membered rings
        
        # Count number of nitrogen atoms in the ring
        n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_nitrogen != 1:
            continue
        
        # Collect all bonds that connect atoms in the ring (intra-ring bonds)
        ring_atom_set = set(ring)
        intra_ring_bonds = []
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_atom_set and a2 in ring_atom_set:
                intra_ring_bonds.append(bond)
        
        if len(intra_ring_bonds) != 5:
            # Incomplete or weird ring – skip it.
            continue
        
        # Count explicit double bonds among intra-ring bonds.
        intra_double = sum(1 for bond in intra_ring_bonds if bond.GetBondType() == rdchem.BondType.DOUBLE)
        
        # Now also count “exocyclic” multiple bonds.
        # For each atom in the ring, examine bonds to atoms not in the ring.
        exo_double = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                # if the neighbor is not in the ring
                nb = bond.GetOtherAtom(atom)
                if nb.GetIdx() in ring_atom_set:
                    continue
                # We count if the bond is a double bond.
                if bond.GetBondType() == rdchem.BondType.DOUBLE:
                    # For our purpose, only count exocyclic multiple bonds to heteroatoms O or S.
                    if nb.GetAtomicNum() in (8, 16):
                        exo_double += 1
        total_multiple = intra_double + exo_double
        
        # For a dihydropyrrole core, we expect exactly one multiple bond.
        if total_multiple == 1:
            return True, "Found a 5-membered ring with one nitrogen and one multiple bond (dihydropyrrole core)."
    
    return False, "No 5-membered pyrroline (dihydropyrrole) core meeting our criteria found."

# Example usage (for testing):
if __name__ == '__main__':
    # Test one known pyrroline: 2-acetyl-1-pyrroline
    test_smiles = "CC(=O)C1=NCCC1"
    result, reason = is_pyrroline(test_smiles)
    print(result, reason)