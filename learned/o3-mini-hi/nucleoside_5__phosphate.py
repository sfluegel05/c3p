"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: Nucleoside 5'-phosphate
Definition: A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base 
in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.
Heuristics used here:
  - At least one nucleobase is detected (an aromatic ring with ≥2 nitrogen atoms).
  - A sugar (ribose/deoxyribose) is detected either as a five-membered ring (4C+1O)
    or as an open-chain fragment showing a CH2 group attached to a phosphate.
  - The sugar part is directly connected to the nucleobase.
  - A CH2 on the sugar (or sugar-like fragment) is attached via oxygen to a phosphorus,
    corresponding to the phosphorylated C-5 position.
  - The overall molecular weight is not excessive (here we require <=800 Da) to avoid
    classifying large nucleotide-containing cofactors (e.g. CoA derivatives).
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
    
    # As a rough check, ensure the molecule is not huge.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 800:
        return False, f"Molecular weight {mol_wt:.1f} too high for a typical nucleoside phosphate (<800 Da)"
        
    # --------- Step 1. Look for an aromatic nucleobase candidate -----------
    # We assume a nucleobase is an aromatic ring containing at least 2 nitrogen atoms.
    ring_info = mol.GetRingInfo()
    nucleobase_atom_indices = set()
    nucleobase_found = False
    for ring in ring_info.AtomRings():
        # Check if all atoms in this ring are aromatic and count nitrogens.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            n_nitrogens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogens >= 2:
                nucleobase_found = True
                nucleobase_atom_indices.update(ring)
                # We consider the first candidate ring found.
                break
    if not nucleobase_found:
        return False, "No nucleobase candidate (aromatic ring with ≥2 nitrogens) found"
    
    # ---------- Step 2. Find the sugar (ribose/deoxyribose) candidate -----------
    sugar_ring_found = False
    sugar_ring_atoms = set()
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            n_oxygens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            n_carbons = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if n_oxygens == 1 and n_carbons == 4:
                sugar_ring_found = True
                sugar_ring_atoms = set(ring)
                break
    sugar_candidate_found = False  # will be True if either closed or open candidate found
    candidate_atoms = None          # indices making up the sugar (ring or chain)
    
    if sugar_ring_found:
        sugar_candidate_found = True
        candidate_atoms = sugar_ring_atoms
    else:
        # Try to find an open-chain sugar-like fragment by searching for a CH2 group 
        # that is attached via an oxygen to a phosphorus.
        # We use a SMARTS pattern for CH2-O-P.
        open_chain_smarts = "[CH2]-[O]-[P]"
        patt = Chem.MolFromSmarts(open_chain_smarts)
        if mol.HasSubstructMatch(patt):
            # For our purposes, we take the CH2 group(s) of the match and then
            # include their immediate neighbors (as the sugar candidate).
            matches = mol.GetSubstructMatches(patt)
            # Choose one match (if many)
            ch2_idx, oxy_idx, p_idx = matches[0]
            candidate_atoms = {ch2_idx, oxy_idx, p_idx}
            sugar_candidate_found = True
    if not sugar_candidate_found:
        return False, "No sugar (closed 5-membered ring or open-chain CH2-O-P fragment) found"

    # --------- Step 3. Check that the sugar candidate is attached to the nucleobase -----------
    # We require that at least one atom from the sugar candidate is directly bonded to
    # an atom from the nucleobase candidate.
    attachment_found = False
    for idx in candidate_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in nucleobase_atom_indices:
                attachment_found = True
                break
        if attachment_found:
            break
    if not attachment_found:
        return False, "Sugar candidate not attached to a nucleobase candidate"

    # --------- Step 4. Check for a phosphate group attached to the sugar (expected at C-5) -----------
    # We look for an sp3 carbon (CH2) within the sugar candidate (if ring) or in its vicinity that has an
    # oxygen neighbor which itself is bound to a phosphorus.
    phosphate_found = False
    # We iterate over atoms in the whole molecule that are carbons (atomic number 6) and check if:
    #   (a) They have at least 2 bound hydrogens
    #   (b) They are adjacent to an atom in our sugar candidate (or candidate extension)
    #   (c) They have a neighbor oxygen that is bound to a phosphorus.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # use GetTotalNumHs() to quickly check for CH2 (at least 2 H)
        if atom.GetTotalNumHs() < 2:
            continue
        # Check if this atom is attached to the sugar candidate
        attached_to_sugar = any(nbr.GetIdx() in candidate_atoms for nbr in atom.GetNeighbors())
        if not attached_to_sugar:
            continue
        # Now look at its neighbors: want an oxygen that bonds to a phosphorus.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                for nn in nbr.GetNeighbors():
                    if nn.GetAtomicNum() == 15:  # phosphorus
                        phosphate_found = True
                        break
                if phosphate_found:
                    break
        if phosphate_found:
            break
    if not phosphate_found:
        return False, "No phosphate group found attached (via oxygen) to a CH2 (expected at C-5) of the sugar fragment"
    
    # --------- All criteria are met -----------
    return True, "Molecule contains a nucleobase (aromatic heterocycle with ≥2 N), a sugar (ribose/deoxyribose or open-chain sugar-like fragment) attached to it, and a phosphate at the likely C-5 position"

# (Optional) Testing examples:
# test_smiles = "O[C@@H]1[C@@H](COP(O)(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O"  # uridine 5'-monophosphate
# result, reason = is_nucleoside_5__phosphate(test_smiles)
# print(result, reason)