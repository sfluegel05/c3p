"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: A cyclitol carboxylic acid (quinic acid and derivatives)

A quinic acid core is defined as:
  - A six-membered ring made entirely of nonaromatic, sp3-hybridized carbons in which all bonds are single.
  - The ring must have exactly five substituents attached (neighbors that are not part of the ring). 
    Each substituent must be one of:
      (a) An oxygen atom (which should represent a hydroxyl or acyloxy link), or 
      (b) A carboxyl group attached via a carbon that shows a double bond to an oxygen (the “free” carboxyl group).
  - In addition, at least one of these substituents must unambiguously be a free carboxyl.
  
This approach attempts to be more selective to avoid false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid derivative based on its SMILES string.
    
    The criteria are:
      1. It must contain at least one six-membered ring in which every atom is a nonaromatic sp3 carbon,
         with all bonds between ring atoms being single.
      2. One of these rings must be substituted with exactly five substituents (neighbors that are not in the ring).
         Each substituent must be either:
            - a single oxygen (which we assume comes from an –OH or an acyloxy group),
            - or a carbon that qualifies as a free carboxyl group.
      3. At least one of the substituents must be a free carboxyl group.
    
    For a carbon substituent to be considered a free carboxyl group, we require that:
         - It has exactly three neighbors (one is the ring carbon) and exactly two oxygens.
         - One of the bonds from the substituent carbon to oxygen is a double bond,
           and the oxygen in the single bond must have at least one hydrogen (indicating a free –COOH,
           rather than an ester carbonyl which often has the oxygen linked further).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a quinic acid derivative, False otherwise.
        str: Reason for the classification.
    """
    # Parse molecule and add explicit hydrogens (to better catch O-H groups)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Define helper function to test if a given atom (a neighbor of the ring) qualifies as a free carboxyl group.
    def is_free_carboxyl(atom):
        # Must be carbon
        if atom.GetAtomicNum() != 6:
            return False
        # Get neighbors of the carbon (one should be the ring; the rest should be oxygens)
        nbrs = atom.GetNeighbors()
        # Expect exactly 3 neighbors: one from the ring and two oxygens.
        if len(nbrs) != 3:
            return False
        oxy_count = 0
        dbl_to_o = False
        single_to_o_with_H = False
        for nbr in nbrs:
            if nbr.GetAtomicNum() == 8:
                oxy_count += 1
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                if bond.GetBondType() == rdchem.BondType.DOUBLE:
                    dbl_to_o = True
                elif bond.GetBondType() == rdchem.BondType.SINGLE:
                    # Check if this oxygen has at least one explicit hydrogen (or is deprotonated,
                    # i.e. has a negative formal charge) – likely indicating a free carboxyl.
                    for o_nbr in nbr.GetNeighbors():
                        if o_nbr.GetAtomicNum() == 1:
                            single_to_o_with_H = True
                            break
                    # Alternatively, if the oxygen carries a negative charge it may be deprotonated.
                    if nbr.GetFormalCharge() < 0:
                        single_to_o_with_H = True
            else:
                # If a neighbor is not oxygen and is not the ring atom, disqualify.
                # (We expect exactly one non-oxygen neighbor, which is the ring atom.)
                pass
        # Should have exactly 2 oxygens, one double-bond and one single-bond with hydrogen.
        if oxy_count == 2 and dbl_to_o and single_to_o_with_H:
            return True
        return False

    # Loop through all rings and check if any qualifies as a quinic acid core.
    for ring in rings:
        # Only consider six-membered rings.
        if len(ring) != 6:
            continue
        
        # Check that each atom in the ring is a nonaromatic sp3 hybridized carbon.
        valid_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                valid_ring = False
                break
            if atom.GetIsAromatic():
                valid_ring = False
                break
            if atom.GetHybridization() != rdchem.HybridizationType.SP3:
                valid_ring = False
                break
        if not valid_ring:
            continue

        # Ensure all bonds between adjacent ring atoms are single bonds.
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

        # Now collect unique substituents attached to the ring.
        substituent_atoms = {}
        # For each atom in the ring, check its neighbors that are not in the ring.
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Use neighbor index as key to avoid double‐counting substituents that connect to two ring atoms.
                substituent_atoms[nbr.GetIdx()] = nbr
                
        if len(substituent_atoms) != 5:
            # The candidate ring must be decorated with exactly 5 substituents.
            continue

        valid_substituent_count = 0
        free_carboxyl_found = False

        # Evaluate each substituent.
        for sub in substituent_atoms.values():
            # Case 1: If the substituent is oxygen then assume it is from an –OH or acyloxy linkage.
            if sub.GetAtomicNum() == 8:
                valid_substituent_count += 1
            # Case 2: If the substituent is carbon, check if it qualifies as a carboxyl group.
            elif sub.GetAtomicNum() == 6:
                if is_free_carboxyl(sub):
                    valid_substituent_count += 1
                    free_carboxyl_found = True
            else:
                # Any other substituent disqualifies this ring.
                pass

        if valid_substituent_count == 5 and free_carboxyl_found:
            msg = ("Found six-membered saturated cyclohexane ring with exactly 5 substituents "
                   "(hydroxy or acyloxy/carboxyl) and at least one free carboxyl group, consistent with a quinic acid core")
            return True, msg

    return False, "No cyclohexane ring with the required five substituents and free carboxyl group detected"

# Example test cases (uncomment these lines to try them):
# test_smiles_list = [
#     "O[C@H]1C[C@@](O)(C[C@H](O)[C@H]1O)C(O)=O",  # (+)-quinic acid (should be True)
#     "C1=CC(=C(C=C1)O)O"  # Hydroquinone (should be False)
# ]
# for smi in test_smiles_list:
#     result, reason = is_quinic_acid(smi)
#     print(smi, result, reason)