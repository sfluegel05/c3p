"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: quinone
Definition: Compounds having a fully conjugated cyclic dione structure, such as that of benzoquinones,
derived from aromatic compounds by conversion of an even number of -CH= groups into -C(=O)- groups,
with any necessary rearrangement of double bonds. (Polycyclic and heterocyclic analogues are included).

Our improved approach:
  - Parses the molecule from a SMILES string.
  - Gets all rings from the molecule.
  - For each candidate ring that is exactly six atoms in size we:
       (a) check that every atom is aromatic or sp2-hybridized,
       (b) check that every bond connecting consecutive atoms in the ring is conjugated,
       (c) identify exocyclic carbonyl groups (a double bond from a ring atom to oxygen not in the ring).
  - If at least two carbonyls are found, and in the case of exactly two groups in a six‐membered ring they are “para‐related” (cyclic distance 3),
    the ring is considered quinone-like.
  - If any such ring is found the molecule is classified as a quinone.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    
    Strategy:
      1. Parse the SMILES string.
      2. Get all ring atom sets from the molecule.
      3. Restrict attention to rings of exactly 6 atoms (the typical size for quinone rings).
      4. For each such ring, check that:
            - All atoms are aromatic or sp2‐hybridized.
            - Every bond between consecutive ring atoms (assuming cyclic order from ring info)
              is conjugated.
            - There are at least two exocyclic carbonyl groups (a C=O on a ring atom with the O out‐of–ring).
            - And if exactly two carbonyl groups are present, that these positions are para‐related
              (i.e. their indices in the ring list differ by 3 when taken cyclically).
      5. If any ring satisfies these conditions, classify the molecule as a quinone.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a quinone; False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings from the molecule (as tuples of atom indices)
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings found in molecule"
    
    # Iterate over rings; restrict to rings of exactly 6 atoms.
    for ring in rings:
        if len(ring) != 6:
            continue  # only consider six–membered rings
        
        # Check that every atom in the ring is aromatic or sp2–hybridized.
        fully_conjugated = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not (atom.GetIsAromatic() or atom.GetHybridization() == rdchem.HybridizationType.SP2):
                fully_conjugated = False
                break
        if not fully_conjugated:
            continue
        
        # Now, check that every bond connecting ring atoms (in cyclic order) is conjugated.
        # Note: The order in ring tuple is assumed to be in-cycle.
        bonds_conjugated = True
        ring_order = list(ring)
        ring_size = len(ring_order)
        for i in range(ring_size):
            j = ring_order[(i + 1) % ring_size]
            bond = mol.GetBondBetweenAtoms(ring_order[i], j)
            if bond is None or not bond.GetIsConjugated():
                bonds_conjugated = False
                break
        if not bonds_conjugated:
            continue
        
        # Identify exocyclic carbonyl groups:
        # For each atom in the ring, if it forms a double bond to an oxygen not in the ring, count it.
        carbonyl_positions = []  # positions (0-indexed within the ring list) with exocyclic C=O
        for pos, idx in enumerate(ring_order):
            atom = mol.GetAtomWithIdx(idx)
            has_ext_carbonyl = False
            for bond in atom.GetBonds():
                # Consider only double bonds as potential carbonyl bonds.
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                neighbor = bond.GetOtherAtom(atom)
                # Look for oxygen neighbor that is NOT a member of the ring
                if neighbor.GetAtomicNum() != 8:
                    continue
                if neighbor.GetIdx() in ring:
                    continue
                has_ext_carbonyl = True
                break
            if has_ext_carbonyl:
                carbonyl_positions.append(pos)
        
        # Must have at least two exocyclic carbonyl groups to be a quinone.
        if len(carbonyl_positions) < 2:
            continue
        
        # For six-membered rings with exactly two exocyclic carbonyls, check that they are para–related.
        if len(carbonyl_positions) == 2:
            pos1, pos2 = carbonyl_positions
            d = abs(pos1 - pos2)
            cyclic_distance = min(d, ring_size - d)
            if cyclic_distance != 3:
                continue
        
        reason = (f"Found quinone ring with atoms {tuple(ring_order)}: "
                  f"{len(carbonyl_positions)} exocyclic carbonyl group(s) on a six-membered fully conjugated ring.")
        return True, reason

    return False, "No six-membered fully conjugated ring with at least two exocyclic carbonyl groups (with proper para–arrangement) was found."

# Example calls (uncomment to test)
# print(is_quinone("Oc1ccc2C(=O)c3c(O)ccc(O)c3C(=O)c2c1O"))  # quinalizarin (expected: True)
# print(is_quinone("O=C1C2=C(O)C=C(C)C=C2C(=O)[C@@]34C15C(=O)[C@H]6C(OC)C3C7[C@]86C5(C(=O)C9=C(O)C=C(C)C=C9C8=O)C(=O)[C@H]4C7OC"))  # Rugulin (expected: False or not quinone)