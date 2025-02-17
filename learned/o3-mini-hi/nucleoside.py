"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: Nucleoside
Definition: An N-glycosyl compound that has both a nucleobase 
           (normally adenine, guanine, xanthine, thymine, cytosine or uracil)
           and either a ribose or deoxyribose as functional parents.
Improvements:
    - Reject molecules with phosphorus (e.g. nucleotides) or with very high molecular weight (>500 Da).
    - Identify a candidate sugar ring: a five-membered ring (furanose) containing exactly 1 oxygen and 4 carbons.
    - Instead of strictly enforcing that a substituent has exactly two hydrogens, we now require that some sugar carbon
      bears an exocyclic neighbor which is a carbon that in turn is attached to at least one oxygen bearing a hydrogen.
    - Identify an N-glycosidic bond between a sugar carbon (typically the anomeric carbon) and a nucleobase nitrogen.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    
    A nucleoside is defined as an N-glycosyl compound that contains:
      - A five-membered sugar ring (ribose or deoxyribose) with exactly one oxygen and four carbons.
      - An exocyclic substituent on the sugar ring corresponding to the 5′-group,
        which is normally a CH2OH (or related) fragment.
      - A nucleobase (typically a 5- or 6-membered heterocycle containing at least one nitrogen)
        that is bonded via an N-glycosidic bond to a sugar carbon.
    
    Additional filters:
      - Molecules containing phosphorus atoms (likely nucleotides) are rejected.
      - Molecules with high molecular weight (>500 Da) are rejected.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): True along with a message if the molecule qualifies as a nucleoside,
                     False along with the reason if it does not.
    """
    # Parse SMILES and check it is valid.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules with phosphorus (likely nucleotides)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus (likely a nucleotide rather than a nucleoside)"
    
    # Reject molecules that are too heavy (nucleosides are usually less than 500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) for a simple nucleoside"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()

    sugar_ring = None
    # Look for a five-membered ring with exactly 1 oxygen and 4 carbons (candidate ribose or deoxyribose)
    for ring in all_rings:
        if len(ring) != 5:
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        num_oxygen = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        num_carbon = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if num_oxygen == 1 and num_carbon == 4:
            sugar_ring = set(ring)
            break
    if sugar_ring is None:
        return False, ("No ribose or deoxyribose sugar ring found " 
                       "(expected a five-membered ring with exactly 1 oxygen and 4 carbons)")
    
    # Add explicit hydrogens to have a better handle on hydrogen counts.
    mol_h = Chem.AddHs(mol)
    
    # Check for an exocyclic substituent on the sugar ring.
    # We search for a neighbor of a sugar carbon (in sugar_ring) that is a carbon (outside the ring)
    # and that carbon in turn is attached to at least one oxygen (also outside the ring) that bears a hydrogen.
    exocyclic_found = False
    for idx in sugar_ring:
        atom = mol_h.GetAtomWithIdx(idx)
        # Only consider carbon atoms of the sugar ring.
        if atom.GetAtomicNum() != 6:
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in sugar_ring:
                continue  # avoid atoms inside the sugar ring
            # Look for a carbon neighbor (candidate for CH2 group)
            if nbr.GetAtomicNum() == 6:
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetIdx() in sugar_ring:
                        continue
                    if subnbr.GetAtomicNum() == 8:
                        # Check that the oxygen has at least one hydrogen attached.
                        if subnbr.GetTotalNumHs() >= 1:
                            exocyclic_found = True
                            break
                if exocyclic_found:
                    break
        if exocyclic_found:
            break
    if not exocyclic_found:
        return False, "Sugar ring missing expected exocyclic CH2OH substituent (5'-group)"
    
    # Identify an N-glycosidic bond between the sugar and the nucleobase.
    # We look for any bond with one end in the sugar ring (and being a carbon)
    # and the other end being a nitrogen that is part of a heterocycle with 5 or 6 atoms.
    glyco_bond_found = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        idx1, idx2 = a1.GetIdx(), a2.GetIdx()
        if (idx1 in sugar_ring and idx2 not in sugar_ring) or (idx2 in sugar_ring and idx1 not in sugar_ring):
            # Identify the sugar atom and the candidate nucleobase atom.
            sugar_atom = a1 if idx1 in sugar_ring else a2
            candidate = a2 if idx1 in sugar_ring else a1
            # The sugar atom should be carbon.
            if sugar_atom.GetAtomicNum() != 6:
                continue
            # The bonded atom should be nitrogen.
            if candidate.GetAtomicNum() != 7:
                continue
            # Let us check that the candidate nitrogen belongs to a ring (of size 5 or 6)
            # that is not wholly part of the sugar ring.
            candidate_rings = []
            for ring in all_rings:
                if candidate.GetIdx() in ring and not set(ring).issubset(sugar_ring):
                    candidate_rings.append(ring)
            if not candidate_rings:
                continue
            # Among these candidate rings, at least one should be a heterocycle (5 or 6 atoms) containing nitrogen.
            valid_base = False
            for ring in candidate_rings:
                if len(ring) in (5, 6):
                    n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
                    if n_count >= 1:
                        valid_base = True
                        break
            if valid_base:
                glyco_bond_found = True
                break
    if not glyco_bond_found:
        return False, "No N-glycosidic bond between sugar and nucleobase detected"
    
    return True, ("Molecule is a nucleoside: contains a five-membered sugar ring (ribose/deoxyribose) with an exocyclic CH2OH (or related) "
                  "group and a nucleobase (5- or 6-membered heterocycle with N) linked via an N-glycosidic bond")

# Example usage:
if __name__ == "__main__":
    # Test with one example SMILES: 1-(beta-D-Ribofuranosyl)-1,2-dihydropyrimidine.
    test_smiles = "O1[C@@H](N2CN=CC=C2)[C@H](O)[C@H](O)[C@H]1CO"
    result, reason = is_nucleoside(test_smiles)
    print(result, ":", reason)