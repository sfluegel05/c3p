"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: Nucleoside 5'-phosphate
Definition: A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base 
in which C-5 of the (deoxy)ribose moiety is phosphorylated (mono-, di-, tri- or tetra-phosphate).
Heuristics used here:
  - The molecule must have moderate size (MW <= 800 Da).
  - A nucleobase is detected as an aromatic ring with at least 2 nitrogen atoms.
  - A sugar candidate is firstly detected as a 5-membered ring (furanose) with 4 carbons and 1 oxygen 
    that is directly connected to the nucleobase (via its anomeric carbon).
  - The sugar candidate must have an exocyclic CH2 group (i.e. the candidate C-5) which is bound (via oxygen)
    to a phosphorus atom.
  - If no closed sugar is found, a fallback using an open-chain pattern [CH2]-O-P is applied, again requiring a nearby nucleobase.
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
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Reject molecules that are very heavy (to avoid large cofactors)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 800:
        return False, f"Molecular weight {mol_wt:.1f} too high for a typical nucleoside phosphate (<800 Da)"
    
    # ------------- Step 1: Identify nucleobase candidate -------------
    # A nucleobase candidate is defined as any aromatic ring with at least 2 nitrogen atoms.
    nucleobase_atoms = set()
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        # Only consider rings where every atom is aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            n_nitrogens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogens >= 2:
                nucleobase_atoms.update(ring)
    if not nucleobase_atoms:
        return False, "No nucleobase candidate (aromatic ring with ≥2 nitrogens) found"
    
    # ------------- Step 2: Identify sugar candidate from a closed (furanose) ring -------------
    # Look for 5-membered rings with exactly 4 carbons and 1 oxygen.
    sugar_candidate = None  # Will be a set of atom indices that form the ring plus possibly a connected exocyclic CH2
    anomeric_atom = None  # atom in the sugar ring attached to nucleobase
    sugar_found = False
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        if len(ring) == 5:
            n_oxygens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            n_carbons = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if n_oxygens == 1 and n_carbons == 4:
                # Check that at least one ring atom is directly connected to a nucleobase atom.
                attached_to_base = False
                candidate_anomeric = None
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() in nucleobase_atoms:
                            attached_to_base = True
                            candidate_anomeric = idx
                            break
                    if attached_to_base:
                        break
                if not attached_to_base:
                    continue  # skip this ring candidate
                
                # Now search for an exocyclic CH2 group that is the C-5 candidate.
                # In a typical nucleoside the sugar ring is attached (by the anomeric carbon)
                # to the nucleobase and one of the OTHER ring carbons should have a neighbor that is a CH2.
                exocyclic_found = False
                for idx in ring:
                    if idx == candidate_anomeric:
                        continue  # skip anomeric carbon
                    atom = mol.GetAtomWithIdx(idx)
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() in ring:
                            continue  # inside ring already
                        # Look for a carbon which might be CH2.
                        if nbr.GetAtomicNum() == 6 and nbr.GetHybridization() in (Chem.rdchem.HybridizationType.SP3,):
                            # Check that it has two hydrogens
                            if nbr.GetTotalNumHs() == 2:
                                # Now check this CH2 (candidate C-5) for a phosphate attachment.
                                phosphate_attached = False
                                for subnbr in nbr.GetNeighbors():
                                    # We want an oxygen directly attached to the CH2 that is bound to phosphorus.
                                    if subnbr.GetAtomicNum() == 8:
                                        # Check if this oxygen has a neighbor with phosphorus.
                                        for oxnbr in subnbr.GetNeighbors():
                                            if oxnbr.GetAtomicNum() == 15:
                                                phosphate_attached = True
                                                break
                                    if phosphate_attached:
                                        break
                                if phosphate_attached:
                                    # We found a sugar ring candidate with a CH2-O-P branch.
                                    sugar_candidate = set(ring)
                                    anomeric_atom = candidate_anomeric
                                    exocyclic_found = True
                                    break
                    if exocyclic_found:
                        break
                if exocyclic_found:
                    sugar_found = True
                    break
    
    # ------------- Step 3: Fallback to open-chain sugar detection -------------
    # In some cases the sugar may be in an open form.
    if not sugar_found:
        # Look for a minimal open-chain fragment: a CH2 group bound via oxygen to phosphorus.
        # We also require that the CH2 has at least one neighbor (other than the oxygen) that is in a nucleobase ring.
        open_chain_smarts = "[CH2]-[O]-[P]"
        patt = Chem.MolFromSmarts(open_chain_smarts)
        matches = mol.GetSubstructMatches(patt)
        fallback_found = False
        for match in matches:
            # match is a tuple of indices: (CH2, O, P)
            ch2_idx = match[0]
            o_idx = match[1]
            # Check that the CH2 has another neighbor that belongs to a nucleobase candidate.
            ch2_atom = mol.GetAtomWithIdx(ch2_idx)
            neighbor_to_base = False
            for nbr in ch2_atom.GetNeighbors():
                if nbr.GetIdx() in nucleobase_atoms:
                    neighbor_to_base = True
                    break
            if neighbor_to_base:
                sugar_found = True
                fallback_found = True
                break
        if not fallback_found:
            return False, "No sugar candidate found: neither a closed 5-membered sugar ring with an exocyclic CH2-O-P branch nor an open-chain CH2-O-P near a nucleobase"

    # ------------- Step 4: Reassurance check: that the sugar candidate is connected to the nucleobase -------------
    # (Already partly done in steps above, but add a final check.)
    attachment_found = False
    if sugar_candidate:
        for idx in sugar_candidate:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in nucleobase_atoms:
                    attachment_found = True
                    break
            if attachment_found:
                break
    if not attachment_found:
        return False, "Sugar candidate is not attached to a nucleobase candidate"

    # ------------- Step 5: Final check: Ensure at least one phosphate group exists attached via oxygen to the sugar candidate -------------
    phosphate_found = False
    # We check all atoms in the molecule that are phosphorus (atomic num 15) and see if they are attached via an O to the sugar candidate or its exocyclic CH2.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8: # oxygen
                    # Check if that oxygen is bound to any atom in the sugar candidate or, in the fallback, to the open-chain CH2.
                    for subnbr in nbr.GetNeighbors():
                        if sugar_candidate and subnbr.GetIdx() in sugar_candidate:
                            phosphate_found = True
                            break
                    if phosphate_found:
                        break
            if phosphate_found:
                break
    if not phosphate_found:
        return False, "No phosphate group found attached (via oxygen) to the sugar candidate"

    return True, ("Molecule contains a nucleobase (aromatic heterocycle with ≥2 N), a sugar (ribose/deoxyribose or sugar-like fragment) "
                  "that is attached to the nucleobase, and a phosphate group linked via an exocyclic CH2 (C-5) of the sugar")
    
# (Optional) Examples for testing -- Uncomment to test:
# test_smiles_list = [
#     "CC(=O)Nc1ccn([C@@H]2O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]2O)c(=O)n1",  # N(4)-acetylcytidine 5'-monophosphate
#     "C[C@@H](OP(O)(=O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O)",  # FMN-L-threonine
#     "NC1=NC(=O)N(C=C1)[C@H]1C[C@H](O)[C@@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)O1"  # dCTP
# ]
# for smi in test_smiles_list:
#     result, reason = is_nucleoside_5__phosphate(smi)
#     print(smi, result, reason)