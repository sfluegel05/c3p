"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: triterpenoid saponin
Definition: A terpene glycoside in which the terpene moiety is a triterpenoid.
That is, the molecule must have at least one sugar (glycoside) unit and an aglycone part with roughly 30 carbon atoms
and several fused (non‐aromatic) rings.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops
from rdkit.Chem import SanitizeFlags

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    
    Heuristic:
      1. The SMILES must be valid.
      2. At least one sugar moiety is present. We identify sugar rings as 5- or 6-membered rings 
         that contain exactly one oxygen and the remaining atoms as carbons.
      3. Remove all atoms belonging to a sugar ring. The largest remaining fragment is taken as the aglycone.
      4. The aglycone should have between 25 and 40 carbon atoms (typical for triterpenoids).
      5. The aglycone must contain at least 4 rings and none of these rings should be aromatic.
      6. Overall, the molecule should be heavy (molecular weight > 500 Da).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a triterpenoid saponin, False otherwise.
        str: Reason for classification.
    """
    # 1. Try parsing the SMILES.
    try:
        # First, try with default sanitization.
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string."
    except Exception as e:
        return False, f"Error parsing SMILES: {e}"
    
    # Sometimes molecules fail during sanitization because of kekulization.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        try:
            # Retry sanitization while skipping the kekulization step.
            Chem.SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_KEKULIZE)
        except Exception as ee:
            return False, f"Sanitization error even with kekulization skipped: {ee}"
    
    # 2. Identify sugar rings: rings with 5 or 6 atoms containing exactly one oxygen and the rest carbons.
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
                    others += 1  # any other hetero atoms
            # A typical sugar (pyranose or furanose) has one oxygen and (ring size - 1) carbons.
            if oxy_count == 1 and carbon_count == (len(ring) - 1) and others == 0:
                sugar_rings.append(set(ring))
    
    if not sugar_rings:
        return False, "No sugar (glycoside) moiety detected."
    
    # 3. Remove sugar ring atoms to calculate the aglycone.
    sugar_atom_indices = set()
    for sring in sugar_rings:
        sugar_atom_indices.update(sring)
        
    editable = Chem.EditableMol(mol)
    for idx in sorted(sugar_atom_indices, reverse=True):
        editable.RemoveAtom(idx)
    aglycone_mol = editable.GetMol()
    
    # Sometimes removal creates disconnected fragments; choose the largest fragment.
    fragments = Chem.GetMolFrags(aglycone_mol, asMols=True, sanitizeFrags=True)
    if not fragments:
        return False, "Sugar removal left no aglycone fragment."
    aglycone = max(fragments, key=lambda m: m.GetNumAtoms())
    
    # 4. Count carbon atoms in the aglycone.
    aglycone_carbons = sum(1 for atom in aglycone.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (25 <= aglycone_carbons <= 40):
        return False, f"Aglycone carbon count not in acceptable range (found {aglycone_carbons})."
    
    # 5. Analyze the ring system of the aglycone.
    ring_info_aglycone = aglycone.GetRingInfo()
    aglycone_rings = ring_info_aglycone.AtomRings()
    if len(aglycone_rings) < 4:
        return False, f"Insufficient number of rings in aglycone (found {len(aglycone_rings)})."
    for ring in aglycone_rings:
        if all(aglycone.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            return False, "Aglycone contains an aromatic ring, which is not typical of triterpenoids."
    
    # 6. Overall molecular weight check.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a triterpenoid saponin."
    
    return True, ("Molecule contains a sugar moiety plus a triterpenoid-like aglycone " +
                  f"(aglycone carbons: {aglycone_carbons}, rings: {len(aglycone_rings)}).")

# For testing purposes (these can be removed in production):
if __name__ == "__main__":
    test_smiles = [
        # Known triterpenoid saponin examples
        "OC1(C2C=3C(C4(C(C5(C(CC4)C(C(O)CC5)(C)C)C)CC3)C)(CCC2(CCC1C)C(OC6OC(C(O)C(O)C6O)CO)=O)C",
        "CN(C)c1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](NC(=O)[C@@H]([NH3+])Cc2ccc(O)cc2)[C@H]1O",  # Causes kekulization issues normally
    ]
    
    for smi in test_smiles:
        verdict, reason = is_triterpenoid_saponin(smi)
        print("SMILES:", smi)
        print("Verdict:", verdict, "| Reason:", reason)
        print("-" * 80)