"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: Nucleoside 5'-phosphate
Definition: A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base in which 
the sugar’s C-5 (the exocyclic carbon, typically CH2OH) is phosphorylated (mono-, di-, tri- or tetra-phosphate).

Heuristics used here:
  1. The molecule must be of moderate size (MW <= 800 Da).
  2. A nucleobase candidate is defined as a fully aromatic ring containing at least 2 nitrogen atoms.
  3. A sugar candidate is first searched for as a 5-membered ring (furanose) containing exactly 1 oxygen and 4 carbons.
     • Such a ring must be directly attached to the nucleobase via one of its atoms (the anomeric carbon).
     • In addition, one of the other ring carbons must have an exocyclic substituent: an oxygen that is directly bonded to a phosphorus atom.
     • To avoid confusion (e.g. dinucleotides) the sugar ring should have exactly one connection to the nucleobase.
  4. If no closed sugar ring is found, a fallback looks for a minimal open-chain fragment of type [#6]-[O]-[P] where the carbon is also connected to a nucleobase.
    
If these criteria are met, the function returns True along with a reason explaining the classification.
Otherwise, a False with the reason of failure is returned.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the criteria for nucleoside 5'-phosphate are met, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Reject molecules that are too heavy (likely not simple nucleosides)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 800:
        return False, f"Molecular weight {mol_wt:.1f} is too high for a typical nucleoside phosphate (<800 Da)."
    
    # Step 1: Identify nucleobase candidate.
    # A nucleobase candidate is defined as any fully aromatic ring (all atoms aromatic)
    # with at least 2 nitrogen atoms.
    nucleobase_atoms = set()
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        # Only consider rings where every atom is aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            n_nitrogens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_nitrogens >= 2:
                nucleobase_atoms.update(ring)
    if not nucleobase_atoms:
        return False, "No nucleobase candidate (aromatic ring with ≥2 nitrogen atoms) found."
    
    # Step 2: Look for a closed sugar candidate.
    sugar_candidate = None  # Set of indices in the sugar ring (if found)
    sugar_has_phosphate = False  # Flag to indicate found exocyclic O-P branch
    base_connections = set()  # Which atoms in the ring connect to nucleobase
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        if len(ring) == 5:
            # Count atoms: we want 4 carbons and 1 oxygen
            n_oxygens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            n_carbons = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if n_oxygens != 1 or n_carbons != 4:
                continue

            # Find connections from this ring to the nucleobase candidate.
            ring_base_conn = set()
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in nucleobase_atoms:
                        ring_base_conn.add(idx)
            # We require exactly one connection (anomeric carbon) to the nucleobase.
            if len(ring_base_conn) != 1:
                continue

            # Now check for an exocyclic phosphorylated group.
            # We exclude the anomeric carbon (the one connecting to the nucleobase).
            anomeric = list(ring_base_conn)[0]
            found_exocyclic = False
            for idx in ring:
                if idx == anomeric:
                    continue
                atom = mol.GetAtomWithIdx(idx)
                # Look at all neighbors not in the ring.
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in ring:
                        continue
                    # If the neighbor is oxygen, check if that oxygen is attached to phosphorus.
                    if nbr.GetAtomicNum() == 8:
                        for oxy_nbr in nbr.GetNeighbors():
                            if oxy_nbr.GetAtomicNum() == 15:
                                found_exocyclic = True
                                break
                    if found_exocyclic:
                        break
                if found_exocyclic:
                    break
            
            if found_exocyclic:
                sugar_candidate = set(ring)
                sugar_has_phosphate = True
                # Use the ring we found; no need to search further.
                break

    # Step 3: Fallback: if no proper closed sugar was found, search for an open-chain fragment.
    fallback_found = False
    if not sugar_candidate:
        # Searching for a fragment of type: carbon - oxygen - phosphorus
        # We use a SMARTS to detect a [#6]-[O]-[P] fragment.
        open_chain_smarts = "[#6]-[O]-[P]"
        patt = Chem.MolFromSmarts(open_chain_smarts)
        matches = mol.GetSubstructMatches(patt)
        for match in matches:
            # match is a tuple: (carbon, oxygen, phosphorus)
            c_idx, o_idx, p_idx = match
            # Check if the carbon (candidate sugar remnant) is attached also to a nucleobase candidate.
            c_atom = mol.GetAtomWithIdx(c_idx)
            attached_to_base = False
            for nbr in c_atom.GetNeighbors():
                if nbr.GetIdx() in nucleobase_atoms:
                    attached_to_base = True
                    break
            if attached_to_base:
                sugar_candidate = {c_idx}
                sugar_has_phosphate = True  # by definition of the substructure match
                fallback_found = True
                break

    if not sugar_candidate or not sugar_has_phosphate:
        return False, ("No sugar candidate found that is attached to a nucleobase and bears an exocyclic O-P branch "
                       "(i.e. phosphate on the probable C-5 position).")
    
    # Step 4: Additional reassurance that the sugar candidate is attached to the nucleobase.
    # In the closed ring case, we already ensured that exactly one ring atom connects to the nucleobase.
    # In the fallback, check that the candidate carbon has at least one nucleobase neighbor.
    if not any(nbr.GetIdx() in nucleobase_atoms for idx in sugar_candidate for nbr in mol.GetAtomWithIdx(idx).GetNeighbors()):
        return False, "Sugar candidate is not attached to a nucleobase."
    
    return True, ("Molecule contains a nucleobase (aromatic heterocycle with ≥2 nitrogens) attached to a sugar "
                  "candidate (a 5-membered furanose or open-chain fragment) that has an exocyclic oxygen bound to phosphorus.")
    
# (Optional) Examples for testing -- Uncomment the lines below to test:
# test_smiles = [
#     "CC(=O)Nc1ccn([C@@H]2O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]2O)c(=O)n1",  # N(4)-acetylcytidine 5'-monophosphate (should be True)
#     "C[C@@H](OP(O)(=O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O)",  # FMN-L-threonine (should be True)
#     "NC1=NC(=O)N(C=C1)[C@H]1C[C@H](O)[C@@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)O1",  # dCTP (should be True or borderline)
#     "NCCCCCCNc1ncnc2n(cnc12)[C@@H]1O[C@@H]2COP(O)(=O)O[C@H]2[C@H]1O",  # A dinucleotide-like or extra phosphate derivative (should be False)
# ]
# for smi in test_smiles:
#     result, reason = is_nucleoside_5__phosphate(smi)
#     print(smi, "==>", result, reason)