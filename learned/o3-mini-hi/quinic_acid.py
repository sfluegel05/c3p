"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: A cyclitol carboxylic acid (quinic acid and derivatives)

We require that a molecule have a quinic acid–like core. For our purposes this is defined as:
  – A saturated, non‐aromatic cyclohexane ring: six sp³ carbon atoms with only single bonds between them.
  – At least one carboxyl (or esterified carboxyl) group attached directly to the ring.
  – A total of between 3 and 5 oxygen–containing substituents attached to the ring.
  
Note: We add hydrogens so that hydroxyl groups are explicitly visible.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid derivative based on its SMILES string.
    
    For our purposes quinic acid (and derivatives) must have a six‐membered ring that is:
      - Fully saturated (all sp3, no aromatic atoms or double bonds within the ring).
      - Carbons only.
      - Attacked by oxygen–substituents: at least one of these substituents must be a carboxyl group 
        (free acid or esterified) and the total count of oxygen substituents (typically hydroxy or acylated –OH)
        attached to the ring (non‐ring neighbors only) should be between 3 and 5.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a quinic acid derivative, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES and add explicit hydrogens for easier detection of -OH groups.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Get ring information (each ring is returned as a tuple of atom indices)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Loop over candidate rings
    for ring in rings:
        # Only consider six-membered rings
        if len(ring) != 6:
            continue
        
        # Check that every atom in the candidate ring is a carbon atom, is not aromatic,
        # and that its hybridization is sp3.
        candidate_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                candidate_ring = False
                break
            if atom.GetIsAromatic():
                candidate_ring = False
                break
            # Make sure its hybridization is sp3 (if not, it might participate in a double bond)
            if atom.GetHybridization() != rdchem.HybridizationType.SP3:
                candidate_ring = False
                break
        if not candidate_ring:
            continue

        # Additionally ensure that all bonds between ring atoms are single bonds.
        bonds_ok = True
        # The ring is a set; check every pair of atoms that are connected in the ring graph.
        for i in range(len(ring)):
            a1 = ring[i]
            # v is the next index (wrap-around)
            a2 = ring[(i+1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
                bonds_ok = False
                break
        if not bonds_ok:
            continue

        # For this candidate ring, count oxygen substituents only on atoms not in the ring.
        oxy_count = 0
        carboxyl_found = False
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            # Look at each neighbor not in the ring:
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # skip atoms inside the ring

                # If neighbor is oxygen, count it as an oxygen substituent.
                if nbr.GetAtomicNum() == 8:
                    # If O is part of a hydroxyl group (has at least one H attached) then count it.
                    # In many ester derivatives the O may not have an explicit H; nevertheless we count it,
                    # but we check later that at least one carboxyl group is present.
                    oxy_count += 1
                # In some cases the carboxyl group is attached via a carbon.
                elif nbr.GetAtomicNum() == 6:
                    # Check if this carbon is likely a carboxyl carbon: it should be connected to at least one oxygen
                    # via a double bond.
                    has_double_bonded_O = False
                    for bond in nbr.GetBonds():
                        # Look for a bond where the other atom is oxygen and the bond is a double bond.
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                            has_double_bonded_O = True
                            break
                    if has_double_bonded_O:
                        carboxyl_found = True
                        oxy_count += 1  # count the carboxyl group as one oxygen substituent.
        
        # We want a minimum total oxygen substituent count in (roughly) the range typical for quinic acid:
        if carboxyl_found and (3 <= oxy_count <= 5):
            return True, ("Found six-membered saturated cyclohexane ring with %d oxygen substituents " 
                          "including a carboxyl/ester group (quinic acid core)" % oxy_count)

    return False, "No cyclohexane ring with required oxygen substituents and carboxyl group detected"

# For debugging or testing, you can un-comment the lines below:
# test_smiles = "O[C@H]1C[C@@](O)(C[C@H](O)[C@H]1O)C(O)=O"  # (+)-quinic acid
# result, reason = is_quinic_acid(test_smiles)
# print(result, reason)