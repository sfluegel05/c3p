"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: A cyclitol carboxylic acid (quinic acid and derivatives)

A quinic acid core is defined as:
  - A six-membered ring in which every atom is a nonaromatic, sp3-hybridized carbon,
    and all bonds between ring atoms are single bonds.
  - The ring must have exactly five substituents attached (neighbors not in the ring)
    that are bound directly via a heteroatom. These substituents should be either:
      (a) An oxygen atom (as in a hydroxyl or an acyloxy group) OR
      (b) A carboxyl group attached directly as a carbon (i.e. a –C(=O)O– unit).
  - At least one of the substituents must be a carboxyl group.
  
This refined approach is designed to avoid false positives where other compounds have a cyclohexane
ring decorated with a few oxygen-containing groups.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid derivative based on its SMILES string.
    
    For our purposes, quinic acid (and many derivatives) must have a six-membered cyclohexane
    ring that is:
      - Fully saturated: all atoms in the ring are sp3-hybridized carbons (nonaromatic) and
        every bond between ring atoms is a single bond.
      - Substituted by exactly five groups attached (only counting neighbors that are not in the ring)
        where each group is either an oxygen atom (from a hydroxyl or acyloxy linkage) OR a carboxyl group.
        In the latter case the ring carbon is directly attached to a carbon that shows a double bond to an oxygen.
      - Containing at least one carboxyl (or esterified carboxyl) substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a quinic acid derivative, False otherwise.
        str: Reason for the classification.
    """
    # Parse the molecule and add explicit hydrogens (to make hydroxyl groups explicit)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)

    # Retrieve ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Loop over rings to select a six-membered candidate ring that fits our criteria.
    for ring in rings:
        if len(ring) != 6:
            continue  # Only consider six-membered rings.

        # Check that every atom in the ring is a nonaromatic, sp3 carbon.
        candidate_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Must be carbon.
            if atom.GetAtomicNum() != 6:
                candidate_ring = False
                break
            # Must not be aromatic.
            if atom.GetIsAromatic():
                candidate_ring = False
                break
            # Must be sp3.
            if atom.GetHybridization() != rdchem.HybridizationType.SP3:
                candidate_ring = False
                break
        if not candidate_ring:
            continue

        # Ensure that all bonds between adjacent ring atoms are single bonds.
        bonds_ok = True
        for i in range(len(ring)):
            a1 = ring[i]
            a2 = ring[(i+1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
                bonds_ok = False
                break
        if not bonds_ok:
            continue

        # For atoms in the ring, count substituents that are attached outside the ring.
        substituent_count = 0
        carboxyl_found = False
        # Loop over each ring atom and its neighbors not in the ring.
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            for nbr in ring_atom.GetNeighbors():
                # Only consider substituents not in the ring.
                if nbr.GetIdx() in ring:
                    continue

                # Case 1: Direct oxygen attachment (e.g. hydroxyl or acyloxy link).
                if nbr.GetAtomicNum() == 8:
                    substituent_count += 1
                    # We do not mark it as a carboxyl here because a carboxyl group is attached via carbon.
                
                # Case 2: Possible carboxyl group attached directly as a carbon.
                elif nbr.GetAtomicNum() == 6:
                    # Check if this carbon is a carboxyl carbon.
                    # We require that it has at least one double bond to an oxygen.
                    has_double_bonded_oxygen = False
                    for bond in nbr.GetBonds():
                        other = bond.GetOtherAtom(nbr)
                        # Look for double bond to oxygen.
                        if other.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                            has_double_bonded_oxygen = True
                            break
                    if has_double_bonded_oxygen:
                        substituent_count += 1
                        carboxyl_found = True

        # For a proper quinic acid core, we require exactly five substituents and at least one carboxyl.
        if substituent_count == 5 and carboxyl_found:
            msg = ("Found six-membered saturated cyclohexane ring with 5 substituents (hydroxy or acyloxy/carboxyl) " +
                   "including at least one carboxyl group (quinic acid core)")
            return True, msg

    return False, "No cyclohexane ring with the required five substituents and carboxyl group detected"

# For testing, you may uncomment the following lines:
# test_smiles = "O[C@H]1C[C@@](O)(C[C@H](O)[C@H]1O)C(O)=O"  # (+)-quinic acid
# result, reason = is_quinic_acid(test_smiles)
# print(result, reason)