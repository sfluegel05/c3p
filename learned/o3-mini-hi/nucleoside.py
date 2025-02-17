"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: Nucleoside
Definition: An N-glycosyl compound that has both a nucleobase 
           (normally adenine, guanine, xanthine, thymine, cytosine or uracil)
           and either a ribose or deoxyribose as functional parents.
Improvement notes:
  - Reject molecules with phosphorus (e.g. nucleotides) or very high molecular weight.
  - Locate a candidate sugar ring as a five-membered ring that contains exactly 1 oxygen and 4 carbons.
  - Also require that the sugar ring carries an exocyclic substituent matching a CH2OH (ribose) or CH2O pattern.
  - Identify a bond from a sugar carbon (ideally the anomeric carbon) to a nitrogen atom that is not in the sugar but belongs
    to a separate ring (of size 5 or 6) as a minimal nucleobase.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    
    A nucleoside should include:
      (1) A sugar ring: typically a five-membered (furanose) ring composed of 4 carbons and 1 oxygen.
          In addition, the sugar should carry an exocyclic CH2OH group (the 5' unit).
      (2) A nucleobase: a heterocyclic ring containing at least one nitrogen (typically 5- or 6-membered).
      (3) An N-glycosidic bond: a covalent bond connecting a sugar carbon (usually the anomeric carbon)
          to a nitrogen atom that is part of the nucleobase ring.
    
    Additional filters:
      - Molecules containing phosphorus atoms are rejected (nucleotides or phospho-derivatives are not nucleosides).
      - Molecules with very high molecular weight (>500 Da) are rejected.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule fits as a nucleoside, False otherwise.
        str: Human-readable reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules that contain phosphorus (likely nucleotides or derivatives)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus (likely a nucleotide rather than a nucleoside)"
    
    # Exclude molecules that are unusually heavy (nucleosides are typically <500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight too high ({mol_wt:.1f} Da) for a simple nucleoside"
    
    # Get ring information
    ringInfo = mol.GetRingInfo()
    all_rings = ringInfo.AtomRings()
    
    sugar_ring = None
    # Look for a five-membered ring with exactly one oxygen and four carbons.
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
    
    # Add hydrogens to improve exocyclic group checking later.
    mol_with_h = Chem.AddHs(mol)
    
    # Check that the sugar ring carries an exocyclic CH2OH group.
    # We expect one of the ring carbons to be bonded to an sp3 carbon (CH2) outside the ring,
    # which in turn is bonded to at least one oxygen that itself has at least one hydrogen.
    exocyclic_ok = False
    # Record the index of the sugar atom that has the glycosidic bond later.
    sugar_glyco_idx = None
    for idx in sugar_ring:
        atom = mol_with_h.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:  # we consider only carbon atoms in the ring for exocyclic linkage
            continue
        # Look at neighbors not in the sugar ring.
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in sugar_ring:
                continue
            # We expect the CH2OH group to be sp3 carbon with 2 hydrogens.
            if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                # Check number of attached hydrogens (including implicit)
                if nbr.GetTotalNumHs() == 2:
                    # Now check that this neighbor has at least one oxygen neighbor with a hydrogen.
                    for subnbr in nbr.GetNeighbors():
                        if subnbr.GetAtomicNum() == 8 and subnbr.GetTotalNumHs() >= 1:
                            exocyclic_ok = True
                            break
            if exocyclic_ok:
                break
        if exocyclic_ok:
            break
    if not exocyclic_ok:
        return False, "Sugar ring missing expected exocyclic CH2OH substituent (5'-group)"
    
    # Identify an N-glycosidic bond connecting the sugar and the nucleobase.
    glyco_bond_found = False
    nucleobase_bond_idx = None  # index of the non-sugar atom in the bond
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        idx1, idx2 = a1.GetIdx(), a2.GetIdx()
        # Exactly one of the atoms should be in the sugar ring.
        if (idx1 in sugar_ring and idx2 not in sugar_ring) or (idx2 in sugar_ring and idx1 not in sugar_ring):
            sugar_atom = a1 if idx1 in sugar_ring else a2
            candidate = a2 if idx1 in sugar_ring else a1
            # The sugar atom should be carbon.
            if sugar_atom.GetAtomicNum() != 6:
                continue
            # The candidate for the nucleobase must be nitro­gen.
            if candidate.GetAtomicNum() != 7:
                continue
            # Check that the candidate nitrogen is part of a ring outside the sugar ring.
            cand_rings = []
            for ring in all_rings:
                if candidate.GetIdx() in ring and not set(ring).issubset(sugar_ring):
                    cand_rings.append(ring)
            if not cand_rings:
                continue
            valid_base = False
            for ring in cand_rings:
                if len(ring) in (5, 6):
                    # As a minimal requirement, ensure the ring has at least one nitrogen.
                    n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
                    if n_count >= 1:
                        valid_base = True
                        break
            if not valid_base:
                continue
            glyco_bond_found = True
            sugar_glyco_idx = sugar_atom.GetIdx()
            nucleobase_bond_idx = candidate.GetIdx()
            break
    if not glyco_bond_found:
        return False, "No N-glycosidic bond between sugar and nucleobase detected"
    
    return True, ("Molecule is a nucleoside: contains a five-membered sugar ring with an exocyclic CH2OH group and "
                  "a nucleobase (heterocycle with N) linked via an N-glycosidic bond")

# Example usage:
if __name__ == "__main__":
    # Test with one example SMILES: 1-(beta-D-Ribofuranosyl)-1,2-dihydropyrimidine.
    test_smiles = "O1[C@@H](N2CN=CC=C2)[C@H](O)[C@H](O)[C@H]1CO"
    result, reason = is_nucleoside(test_smiles)
    print(result, ":", reason)