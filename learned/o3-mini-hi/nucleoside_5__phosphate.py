"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: Nucleoside 5'-phosphate
Definition: A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base 
in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.
Heuristics used here:
  - The molecule must be of moderate size (molecular weight <=800 Da).
  - A nucleobase is detected as any aromatic ring that contains at least 2 nitrogen atoms.
    (All candidate rings are aggregated.)
  - A sugar candidate is detected either as:
       (a) a five-membered ring with exactly one oxygen and four carbons and which is directly bonded to the nucleobase, or
       (b) an open‐chain fragment containing a CH2 group that is linked via oxygen to a phosphorus atom.
  - Finally, a phosphate group (P atom attached via oxygen) must be found in the vicinity of the sugar candidate.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule matches the criteria for a nucleoside 5'-phosphate, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules that are too heavy (to avoid large cofactors).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 800:
        return False, f"Molecular weight {mol_wt:.1f} too high for a typical nucleoside phosphate (<800 Da)"
    
    # ---------------- Step 1: Identify the nucleobase candidate ----------------
    # We assume a nucleobase is an aromatic ring with at least 2 nitrogen atoms.
    nucleobase_atoms = set()
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        # Only consider rings in which all atoms are aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            # Count nitrogen atoms in the ring.
            n_nitrogens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogens >= 2:
                nucleobase_atoms.update(ring)
    if not nucleobase_atoms:
        return False, "No nucleobase candidate (aromatic ring with ≥2 nitrogens) found"
    
    # ---------------- Step 2: Locate the sugar candidate ----------------
    sugar_candidate = None
    # First, try to detect a closed (ring) sugar candidate.
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            n_oxygens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            n_carbons = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            # Ribose or deoxyribose ring should have one oxygen and four carbons.
            if n_oxygens == 1 and n_carbons == 4:
                # Check that at least one atom of this ring is directly attached to the nucleobase.
                attached_to_base = False
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() in nucleobase_atoms:
                            attached_to_base = True
                            break
                    if attached_to_base:
                        break
                if attached_to_base:
                    sugar_candidate = set(ring)
                    break
    # If no ring sugar is found, try an open-chain sugar-like candidate.
    if sugar_candidate is None:
        # Look for a CH2 group connected via oxygen to a phosphorus.
        # The SMARTS pattern below targets a CH2 group (-[CH2]-[O]-[P])
        open_chain_smarts = "[CH2]-[O]-[P]"
        patt = Chem.MolFromSmarts(open_chain_smarts)
        if patt is not None and mol.HasSubstructMatch(patt):
            # Pick the first match; note that this is a minimal fragment.
            match = mol.GetSubstructMatches(patt)[0]
            sugar_candidate = set(match)
        else:
            return False, "No sugar candidate (closed 5-membered ring or open-chain CH2-O-P fragment) found"
    
    # Verify again that the sugar candidate is attached to the nucleobase.
    attachment_found = False
    for idx in sugar_candidate:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in nucleobase_atoms:
                attachment_found = True
                break
        if attachment_found:
            break
    if not attachment_found:
        return False, "Sugar candidate not attached to a nucleobase candidate"
    
    # ---------------- Step 3: Verify the presence of a phosphate group ----------------
    # We look for a phosphorus atom – attached via at least one oxygen – in the vicinity of the sugar candidate.
    phosphate_found = False
    # Here we iterate over the atoms in our sugar candidate and check each oxygen neighbor.
    for idx in sugar_candidate:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # Look for oxygen neighbor.
            if nbr.GetAtomicNum() == 8:
                # Check if this oxygen is linked to a phosphorus.
                for nn in nbr.GetNeighbors():
                    if nn.GetAtomicNum() == 15:
                        phosphate_found = True
                        break
                if phosphate_found:
                    break
        if phosphate_found:
            break
    if not phosphate_found:
        return False, "No phosphate group found attached (via oxygen) to the sugar candidate"
    
    # ---------------- All criteria are met ----------------
    return True, ("Molecule contains a nucleobase (aromatic heterocycle with ≥2 N), a sugar (ribose/deoxyribose or sugar-like open-chain fragment) "
                  "directly attached to it, and a phosphate group (likely at the C-5 position)")

# (Optional) Testing examples – Uncomment below to test specific SMILES:
# test_smiles = "O[C@@H]1[C@@H](COP(O)(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O"  # uridine 5'-monophosphate
# result, reason = is_nucleoside_5__phosphate(test_smiles)
# print(result, reason)