"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: triterpenoid saponin
Definition: A terpene glycoside in which the terpene moiety is a triterpenoid.
That is, the molecule must have at least one sugar (glycoside) unit and an aglycone part
with roughly 30 carbon atoms and several fused (non‐aromatic) rings.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.

    The heuristic now is:
      1. The molecule must be a valid structure.
      2. At least one sugar moiety must be present. We identify sugar rings as 5- or 6-membered rings
         that contain exactly one oxygen and the rest carbons.
      3. Remove all atoms that belong to any identified sugar ring.
         Then, the largest remaining fragment is taken as the aglycone.
      4. The aglycone should have a carbon count roughly in the range expected for triterpenoids.
         (We require between 25 and 40 carbons.)
      5. The aglycone part should contain at least 4 rings and none of those rings should be aromatic.
      6. Overall the molecule should be heavy (e.g. molecular weight > 500 Da).

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if it is classified as a triterpenoid saponin, False otherwise.
        str: Reason for classification.
    """
    # 1. Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # 2. Identify sugar rings: look for rings that are 5- or 6-membered and consist of one oxygen and the rest carbons.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    sugar_rings = []
    for ring in all_rings:
        if len(ring) in [5, 6]:
            oxy_count = 0
            carbon_count = 0
            others = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                num = atom.GetAtomicNum()
                if num == 8:
                    oxy_count += 1
                elif num == 6:
                    carbon_count += 1
                else:
                    others += 1
            # Typical sugar ring (pyranose or furanose) should have exactly one oxygen and rest carbons.
            if oxy_count == 1 and others == 0:
                sugar_rings.append(set(ring))
    
    if not sugar_rings:
        return False, "No sugar (glycoside) moiety detected."
    
    # Save indices of all atoms that are in any sugar ring.
    sugar_atom_indices = set()
    for sring in sugar_rings:
        sugar_atom_indices.update(sring)
    
    # 3. To better estimate the aglycone, remove all sugar ring atoms from the molecule.
    # Create an editable mol copy and remove atoms in descending order (to preserve indices).
    editable = Chem.EditableMol(mol)
    # sort indices descending so removal does not mess up mapping
    for idx in sorted(sugar_atom_indices, reverse=True):
        editable.RemoveAtom(idx)
    aglycone_mol = editable.GetMol()
    
    # Sometimes removal produces several disconnected fragments.
    # Choose the largest fragment (by number of atoms) as the putative aglycone.
    fragments = Chem.GetMolFrags(aglycone_mol, asMols=True, sanitizeFrags=True)
    if not fragments:
        return False, "Sugar removal left no aglycone fragment."
    # Choose fragment with maximum number of atoms
    aglycone = max(fragments, key=lambda m: m.GetNumAtoms())
    
    # 4. Count carbon atoms in the aglycone fragment.
    aglycone_carbons = sum(1 for atom in aglycone.GetAtoms() if atom.GetAtomicNum() == 6)
    # We require that the aglycone has roughly 30 carbons;
    # allow a margin – here we require between 25 and 40 carbons.
    if not (25 <= aglycone_carbons <= 40):
        return False, f"Aglycone carbons not in acceptable range (found {aglycone_carbons})."
    
    # 5. Check ring system in the aglycone:
    # a) Count rings in the aglycone.
    ring_info_aglycone = aglycone.GetRingInfo()
    aglycone_rings = ring_info_aglycone.AtomRings()
    if len(aglycone_rings) < 4:
        return False, f"Insufficient number of rings in aglycone (found {len(aglycone_rings)})."
    # b) Ensure that none of the rings are aromatic – triterpenoid cores are typically saturated, not aromatic.
    for ring in aglycone_rings:
        # If all atoms in a ring are aromatic, assume the ring is aromatic.
        if all(aglycone.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            return False, "Aglycone contains an aromatic ring, not typical of triterpenoids."
    
    # 6. Check overall molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a triterpenoid saponin."
    
    # Passed all tests.
    return True, ("Molecule contains a sugar moiety plus a triterpenoid-like aglycone "
                  f"(aglycone carbons: {aglycone_carbons}, rings: {len(aglycone_rings)}).")

# For testing purposes (can be removed or commented out in production):
if __name__ == "__main__":
    # Example SMILES strings (both true and some false cases) could be placed here.
    examples = [
        # A known triterpenoid saponin example:
        "OC1(C2C=3C(C4(C(C5(C(CC4)C(C(O)CC5)(C)C)C)CC3)C)(CCC2(CCC1C)C(OC6OC(C(O)C(O)C6O)CO)=O)C",
        # A false positive example (flavonoid glycoside likely):
        "COc1cc(ccc1O)C(=O)OC[C@H]1O[C@@H](Oc2ccc(cc2OC)[C@@H]2Oc3c(OC)cc(CCCO)cc3[C@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O"
    ]
    for smi in examples:
        verdict, reason = is_triterpenoid_saponin(smi)
        print("SMILES:", smi)
        print("Verdict:", verdict, "| Reason:", reason)
        print("-" * 80)