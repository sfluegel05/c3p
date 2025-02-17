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
    - Look for an exocyclic substituent (expected at the 5'-position) that is a CH2OH (or CH2O…) group.
      Rather than enforcing exactly 2 hydrogens, we now allow a small flexibility (>=1 hydrogen)
      and require that the candidate carbon is sp3 and it is directly bonded to an oxygen that itself bears at least one hydrogen.
    - Identify an N-glycosidic bond between a (typically anomeric) sugar carbon and a nucleobase nitrogen.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    
    A nucleoside should include:
      (1) A sugar ring: typically a five-membered (furanose) ring composed of 4 carbons and 1 oxygen.
          In addition, the sugar should carry an exocyclic substituent corresponding to the 5'-group, 
          usually a CH2OH group (or CH2O in certain depictions).
      (2) A nucleobase: a heterocyclic ring (of size 5 or 6) containing at least one nitrogen.
      (3) An N-glycosidic bond: a covalent bond linking a sugar carbon (commonly the anomeric carbon)
          to a nitrogen atom that belongs to the nucleobase.
    
    Additional filters:
      - Molecules containing phosphorus atoms are rejected (these are likely nucleotides or derivatives).
      - Molecules with high molecular weight (>500 Da) are rejected.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule fits as a nucleoside, False otherwise.
        str: Human-readable reason for classification.
    """
    # Parse SMILES and check that it is valid.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules with phosphorus (likely nucleotides or phospho-derivatives).
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus (likely a nucleotide rather than a nucleoside)"
    
    # Exclude molecules that are unusually heavy (nucleosides are typically less than 500 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) for a simple nucleoside"
    
    # Inspect the ring systems.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    
    sugar_ring = None
    # Look for a five-membered ring with exactly one oxygen and four carbons (candidate ribose/deoxyribose).
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
                       "(expected five-membered ring with exactly 1 oxygen and 4 carbons)")
    
    # Add explicit hydrogens to get more reliable hydrogen counts.
    mol_h = Chem.AddHs(mol)
    
    # Check for an exocyclic substituent on the sugar ring.
    # We look for a neighbor of one of the sugar carbons (atom from sugar_ring)
    # which is an sp3 carbon (candidate CH2 group) attached to an oxygen that is not in the ring.
    exocyclic_found = False
    for idx in sugar_ring:
        atom = mol_h.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:  # look only at carbons in the ring
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in sugar_ring:
                continue  # skip atoms in the sugar ring
            if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                # Relax the hydrogen count requirement: allow at least 1 hydrogen.
                if nbr.GetTotalNumHs() >= 1:
                    # Now check that this candidate carbon has at least one oxygen neighbor (exocyclic oxygen).
                    oxygen_attached = False
                    for subnbr in nbr.GetNeighbors():
                        if subnbr.GetIdx() in sugar_ring:
                            continue
                        if subnbr.GetAtomicNum() == 8 and subnbr.GetTotalNumHs() >= 1:
                            oxygen_attached = True
                            break
                    if oxygen_attached:
                        exocyclic_found = True
                        break
            if exocyclic_found:
                break
        if exocyclic_found:
            break
    if not exocyclic_found:
        return False, "Sugar ring missing expected exocyclic CH2OH substituent (5'-group)"
    
    # Identify an N-glycosidic bond between the sugar and the nucleobase.
    glyco_bond_found = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        idx1, idx2 = a1.GetIdx(), a2.GetIdx()
        # One end must be in the sugar ring and the other not.
        if (idx1 in sugar_ring and idx2 not in sugar_ring) or (idx2 in sugar_ring and idx1 not in sugar_ring):
            sugar_atom = a1 if idx1 in sugar_ring else a2
            candidate = a2 if idx1 in sugar_ring else a1
            # The sugar atom should be carbon.
            if sugar_atom.GetAtomicNum() != 6:
                continue
            # The bond should be to a nitrogen.
            if candidate.GetAtomicNum() != 7:
                continue
            # Check that the candidate nitrogen is part of a ring not entirely inside the sugar.
            candidate_rings = []
            for ring in all_rings:
                if candidate.GetIdx() in ring and not set(ring).issubset(sugar_ring):
                    candidate_rings.append(ring)
            if not candidate_rings:
                continue
            # Check that one of these rings is a heterocycle (of size 5 or 6) and includes at least one nitrogen.
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
    
    return True, ("Molecule is a nucleoside: contains a five-membered sugar ring with an exocyclic CH2OH (or related) group and "
                  "a nucleobase (a 5- or 6-membered heterocycle with at least one N) linked via an N-glycosidic bond")

# Example usage:
if __name__ == "__main__":
    # Test with one example SMILES: 1-(beta-D-Ribofuranosyl)-1,2-dihydropyrimidine.
    test_smiles = "O1[C@@H](N2CN=CC=C2)[C@H](O)[C@H](O)[C@H]1CO"
    result, reason = is_nucleoside(test_smiles)
    print(result, ":", reason)