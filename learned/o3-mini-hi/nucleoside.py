"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: Nucleoside
Definition: An N-glycosyl compound that has both a nucleobase 
           (normally adenine, guanine, xanthine, thymine, cytosine or uracil)
           and either a ribose or deoxyribose as functional parents.

Revisions:
  - Reject molecules with phosphorus or with molecular weight >500 Da.
  - Identify a candidate sugar ring as a five-membered ring with exactly 1 oxygen and 4 carbons.
  - Look for an exocyclic substituent on a sugar carbon that is a CH2 group attached to an –OH,
    by checking that the neighboring carbon has two H's and is bonded to an oxygen bearing at least one hydrogen.
  - Identify an N-glycosidic bond between the sugar and a nucleobase (a nitrogen in a 5- or 6-membered heterocycle).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    
    A nucleoside is defined as an N-glycosyl compound that contains:
      - A five-membered sugar ring (ribose or deoxyribose) with exactly one oxygen and four carbons.
      - An exocyclic substituent on the sugar ring corresponding to the 5′-group,
        expected to be a –CH2OH (or related) fragment.
      - A nucleobase (typically a 5- or 6-membered heterocycle containing one or more nitrogen atoms)
        linked via an N-glycosidic bond to a sugar carbon.
    
    Additional filters:
      - Molecules containing phosphorus (likely nucleotides) are rejected.
      - Molecules with molecular weight >500 Da are rejected.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
       (bool, str): True and an explanation if the molecule is a valid nucleoside;
                    False and a reason if the molecule does not meet the definition.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules containing phosphorus (nucleotides usually)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus (likely a nucleotide rather than a nucleoside)"
    
    # Reject if molecular weight is too high for a simple nucleoside.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) for a simple nucleoside"
    
    # Retrieve ring information and look for a five-membered ring with exactly 1 oxygen and 4 carbons.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    sugar_ring = None
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
        return False, "No ribose or deoxyribose sugar ring found (expected a five-membered ring with 1 O and 4 C)"
    
    # Add explicit hydrogens to better assess the hydrogen counts.
    mol_h = Chem.AddHs(mol)
    
    # Revised exocyclic substituent check:
    # For nucleosides, we expect one of the sugar carbons to have an exocyclic CH2 group that is bonded to an OH.
    exocyclic_found = False
    for idx in sugar_ring:
        atom = mol_h.GetAtomWithIdx(idx)
        # We only care about carbon atoms in the sugar ring.
        if atom.GetAtomicNum() != 6:
            continue
        # Look at neighbors outside the sugar ring.
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in sugar_ring:
                continue
            # The candidate exocyclic substituent should be a carbon.
            if nbr.GetAtomicNum() == 6:
                # Check if this carbon appears as a CH2 (i.e. at least two hydrogens attached).
                if nbr.GetTotalNumHs() < 2:
                    continue
                # Now check that this carbon is attached to an oxygen atom that itself has at least one hydrogen.
                ok_oh = False
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetIdx() in sugar_ring:
                        continue
                    if subnbr.GetAtomicNum() == 8 and subnbr.GetTotalNumHs() >= 1:
                        ok_oh = True
                        break
                if ok_oh:
                    exocyclic_found = True
                    break
        if exocyclic_found:
            break
    if not exocyclic_found:
        return False, "Sugar ring missing expected exocyclic CH2OH substituent (5'-group)"
    
    # Identify an N-glycosidic bond between the sugar and the nucleobase.
    glyco_bond_found = False
    # For each bond, if one end is in the sugar ring (should be a carbon) and the other end is a nitrogen,
    # then check if that nitrogen is part of a heterocycle (ring of size 5 or 6 containing nitrogen).
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        in_ring1 = a1.GetIdx() in sugar_ring
        in_ring2 = a2.GetIdx() in sugar_ring
        if in_ring1 ^ in_ring2:  # xor: one in sugar, one outside.
            sugar_atom = a1 if in_ring1 else a2
            candidate = a2 if in_ring1 else a1
            # The sugar atom should be a carbon.
            if sugar_atom.GetAtomicNum() != 6:
                continue
            # Candidate must be a nitrogen.
            if candidate.GetAtomicNum() != 7:
                continue
            # Check if the candidate nitrogen is part of a ring (of size 5 or 6) that is not completely the sugar ring.
            candidate_rings = []
            for ring in all_rings:
                if candidate.GetIdx() in ring and not set(ring).issubset(sugar_ring):
                    candidate_rings.append(ring)
            if not candidate_rings:
                continue
            # Among these rings, at least one should be of size 5 or 6 and contain at least one nitrogen.
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
    
    return True, ("Molecule is a nucleoside: contains a five-membered sugar ring (ribose/deoxyribose) with an "
                  "exocyclic CH2OH (or related) group and a nucleobase linked via an N-glycosidic bond.")

# Example usage:
if __name__ == "__main__":
    # Test an example SMILES: 1-(beta-D-Ribofuranosyl)-1,2-dihydropyrimidine.
    test_smiles = "O1[C@@H](N2CN=CC=C2)[C@H](O)[C@H](O)[C@H]1CO"
    result, reason = is_nucleoside(test_smiles)
    print(result, ":", reason)