"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: quinone
Definition: Compounds having a fully conjugated cyclic dione structure, such as that of benzoquinones,
derived from aromatic compounds by conversion of an even number of -CH= groups into -C(=O)- groups,
with any necessary rearrangement of double bonds. (Polycyclic and heterocyclic analogues are included).

This improved approach:
  - Parses the molecule from SMILES.
  - Iterates over SSSR rings (of size ≥ 5) and only considers rings that are fully conjugated.
    (We require that each atom is sp2‐hybridized or aromatic and that each bond
     connecting ring atoms is conjugated.)
  - For each ring, we find exocyclic carbonyl groups (an sp2 carbon double‐bonded to oxygen, with the oxygen external to the ring).
  - If at least two such carbonyls are found the ring is a candidate.
    For 6-membered rings with exactly two carbonyls we check that the two positions are “para‐related”
    (cyclic distance 3), though if more than two carbonyls are attached this check is waived.
  - If any ring satisfies these conditions, the molecule is classified as a quinone.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    
    Strategy:
      1. Parse the SMILES.
      2. Get the unique SSSR rings.
      3. For each ring with at least 5 atoms:
            - Check that every atom in the ring is either sp2 hybridized or aromatic.
            - Check that every bond between ring atoms is conjugated.
            - Identify exocyclic carbonyl groups attached to ring atoms.
      4. For a candidate ring:
            - If the ring has fewer than two exocyclic carbonyls, skip.
            - If the ring is 6-membered and exactly two carbonyls are found, require that 
              the two carbonyl-bearing positions are “para‐related” (cyclic distance of 3).
      5. Return True along with a descriptive reason if any ring meets these criteria.
         Otherwise, return False.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as quinone, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the list of rings (each ring is a tuple of atom indices)
    rings = list(Chem.GetSymmSSSR(mol))
    if not rings:
        return False, "No rings found in molecule"
    
    # Iterate over candidate rings (only rings of size at least 5)
    for ring in rings:
        if len(ring) < 5:
            continue
        
        ring_indices = list(ring)
        
        # Verify that every atom in the ring is at least sp2-hybridized or aromatic.
        ring_conjugated = True
        for idx in ring_indices:
            atom = mol.GetAtomWithIdx(idx)
            if not (atom.GetIsAromatic() or atom.GetHybridization() == rdchem.HybridizationType.SP2):
                ring_conjugated = False
                break
        if not ring_conjugated:
            continue
        
        # Verify that every bond connecting atoms in the ring is conjugated.
        bonds_conjugated = True
        ring_bonds = []
        for i, idx in enumerate(ring_indices):
            j = ring_indices[(i+1) % len(ring_indices)]
            bond = mol.GetBondBetweenAtoms(idx,j)
            if bond is None or not bond.GetIsConjugated():
                bonds_conjugated = False
                break
        if not bonds_conjugated:
            continue
        
        # Look for exocyclic carbonyl groups on ring atoms.
        # A carbonyl group is detected when a ring atom (typically sp2 C) is double-bonded to an oxygen
        # where that oxygen is not part of this ring.
        carbonyl_positions = []  # positions (indices in the ring_indices list) that bear an external C=O
        for pos, idx in enumerate(ring_indices):
            atom = mol.GetAtomWithIdx(idx)
            has_ext_carbonyl = False
            for bond in atom.GetBonds():
                # Only consider double bonds
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                neighbor = bond.GetOtherAtom(atom)
                # Check if the neighbor is oxygen
                if neighbor.GetAtomicNum() != 8:
                    continue
                # Exclude if the oxygen is also in the ring.
                if neighbor.GetIdx() in ring_indices:
                    continue
                # Additionally, check that the oxygen has a double-bond environment (i.e. a carbonyl).
                has_ext_carbonyl = True
                break
            if has_ext_carbonyl:
                carbonyl_positions.append(pos)
        
        if len(carbonyl_positions) < 2:
            continue
        
        # For 6-membered rings with exactly two carbonyls, ensure a para relationship.
        if len(ring_indices) == 6 and len(carbonyl_positions) == 2:
            pos1, pos2 = carbonyl_positions
            d = abs(pos1 - pos2)
            cyclic_distance = min(d, 6 - d)
            if cyclic_distance != 3:
                # If not in a para-like arrangement, skip this ring.
                continue
                
        # If we reached here the ring is conjugated and has at least two exocyclic carbonyl groups.
        reason = (f"Found ring with atoms {tuple(ring_indices)}: "
                  f"{len(carbonyl_positions)} exocyclic carbonyl group(s) attached to a fully conjugated ring of size {len(ring_indices)}.")
        return True, reason

    return False, "No fully conjugated ring with at least two exocyclic carbonyl groups (with proper arrangement for 6-membered rings) was found."

# Example calls (uncomment to test)
# print(is_quinone("Oc1ccc2C(=O)c3c(O)ccc(O)c3C(=O)c2c1O"))  # quinalizarin (expected true)
# print(is_quinone("O=C1C(=C[C@@H]2C=C[C@@H]3CC(N[C@]43[C@@H]2C(=O)C5=C(C(O)=C(C)C(=C5C4=O)O)C([C@@H](C=C[C@H]([C@@H]1C)O)C)=O)=O)CC(C)C"))  # Ansaseomycin A (may be false)