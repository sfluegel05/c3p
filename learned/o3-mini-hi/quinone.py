"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: quinone
Definition: Compounds having a fully conjugated cyclic dione structure (e.g. benzoquinones and analogues)
derived from aromatic compounds by conversion of an even number of –CH= groups into –C(=O)– groups.
Our improved approach:
  - Parse the molecule from a SMILES string.
  - Identify rings from the molecule.
  - Consider only rings of exactly 6 atoms that are aromatic and consist entirely of carbon atoms.
  - Check that every bond between consecutive ring atoms is conjugated.
  - For each ring, look at each ring atom for an exocyclic carbonyl group (a double bond from a ring carbon 
    to an oxygen not in the ring).
  - Accept the ring if at least 2 such carbonyl groups are present – and if exactly 2, their positions must be 
    para–related (cyclic distance of 3).
If any such ring is found the molecule is classified as a quinone.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    
    The algorithm:
      1. Parses the SMILES.
      2. Gets all rings from the molecule.
      3. For each ring of exactly 6 atoms:
           - Checks that every atom is aromatic and is a carbon.
           - Verifies that each consecutive bond in the ring is conjugated.
           - Identifies exocyclic carbonyl groups (a double bond from a ring atom to an oxygen outside the ring).
           - Requires at least two such carbonyl groups; if exactly two, they must be para–related (cyclic difference 3).
      4. If any six–membered ring passes these tests, the molecule is classified as a quinone.
      
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        (bool, str): A tuple with True and an explanation if a qualifying ring is found;
                     otherwise, False with the reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings (each ring is a tuple of atom indices)
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings found in molecule"
    
    # Iterate over each ring and restrict to rings of size 6 (typical for benzoquinones)
    for ring in rings:
        if len(ring) != 6:
            continue
        
        # Check that every atom in the ring is carbon and aromatic.
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring):
            continue
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Check that every bond between consecutive ring atoms (using the cyclic order from the ring info) is conjugated.
        ring_order = list(ring)
        ring_size = 6
        bonds_conjugated = True
        for i in range(ring_size):
            a1 = ring_order[i]
            a2 = ring_order[(i + 1) % ring_size]
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond is None or not bond.GetIsConjugated():
                bonds_conjugated = False
                break
        if not bonds_conjugated:
            continue
        
        # Identify exocyclic carbonyl groups:
        # For each ring atom, see if it forms a double bond with an oxygen that is NOT part of the ring.
        carbonyl_positions = []
        for pos, idx in enumerate(ring_order):
            atom = mol.GetAtomWithIdx(idx)
            has_ext_carbonyl = False
            for bond in atom.GetBonds():
                # Only consider double bonds as potential carbonyl bonds.
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                neighbor = bond.GetOtherAtom(atom)
                # Check if neighbor is oxygen and is not in the ring.
                if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in ring:
                    has_ext_carbonyl = True
                    break
            if has_ext_carbonyl:
                carbonyl_positions.append(pos)
        
        if len(carbonyl_positions) < 2:
            continue
        
        # If exactly two exocyclic carbonyl groups, require that they are para–related on the ring (cyclic distance 3).
        if len(carbonyl_positions) == 2:
            pos1, pos2 = carbonyl_positions
            d = abs(pos1 - pos2)
            if min(d, ring_size - d) != 3:
                continue
        
        reason = (f"Found quinone ring: ring atoms {tuple(ring_order)} with "
                  f"{len(carbonyl_positions)} exocyclic carbonyl group(s) on a six-membered aromatic carbon ring.")
        return True, reason

    return False, "No six-membered fully conjugated carbon ring with the required exocyclic carbonyl group(s) found."

# Example calls (uncomment to test)
# print(is_quinone("Oc1ccc2C(=O)c3c(O)ccc(O)c3C(=O)c2c1O"))  # quinalizarin; expected: True
# print(is_quinone("O=C1C2=C(O)C(=C(C)C(=C2C(=O)CC1)O)CC(=O)C"))   # example false positive candidate; expected: False