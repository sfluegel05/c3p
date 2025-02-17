"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: Tertiary amine

A tertiary amine is defined as a compound formally derived from ammonia by replacing
all three hydrogen atoms by hydrocarbyl groups. This implementation improves over a previous
attempt by (1) rejecting candidate nitrogen atoms that appear to be part of amide bonds
(i.e. when the nitrogen is directly attached to a carbon that is double-bonded to oxygen)
and (2) verifying that each substituent attached to the candidate displays predominantly
hydrocarbyl character (roughly, at least 50% of its heavy atoms should be carbon).
For non‐ring candidates we additionally require that removal of the candidate disconnects 
the substituents into three distinct fragments.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tertiary_amine(smiles: str):
    """
    Determines whether the input molecule (given as a SMILES string) qualifies as having
    at least one tertiary amine group (a neutral nitrogen with three heavy substituents that are
    largely hydrocarbyl in composition). 
    
    Returns:
        (bool, str): A tuple containing True (with an explanation) if a tertiary amine is found,
                     and False with a reason otherwise.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens so that bond counts are accurate; necessary for later checks.
    mol_with_H = Chem.AddHs(mol)

    # To track original atom indices reliably after modifications, store each atom's original index.
    for atom in mol_with_H.GetAtoms():
        atom.SetProp("orig_idx", str(atom.GetIdx()))
    
    candidates = []
    # Look for candidate nitrogen atoms.
    for atom in mol_with_H.GetAtoms():
        if atom.GetAtomicNum() != 7:  # not nitrogen
            continue
        if atom.GetFormalCharge() != 0:
            continue
        # Exclude aromatic nitrogen atoms (their lone pairs are delocalized)
        if atom.GetIsAromatic():
            continue
        # Must have no hydrogen atoms explicitly attached.
        if any(neigh.GetAtomicNum() == 1 for neigh in atom.GetNeighbors()):
            continue
        # Must have exactly three heavy (non-hydrogen) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 3:
            continue
        candidates.append(atom)

    if not candidates:
        return False, "No candidate tertiary amine group found in the molecule."

    # helper: check if candidate nitrogen is adjacent to a carbonyl (C=O) group.
    def is_adjacent_to_carbonyl(nitrogen):
        for nbr in nitrogen.GetNeighbors():
            # Check only heavy neighbor atoms (and mostly carbons are expected)
            if nbr.GetAtomicNum() == 6:  # carbon
                for bond in nbr.GetBonds():
                    # If this carbon atom is double-bonded to an oxygen (C=O)
                    other = bond.GetOtherAtom(nbr)
                    if other.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 1:
                        return True
        return False

    # helper: given a collection of atoms, return the ratio of carbon atoms among heavy atoms
    def carbon_fraction(atom_list):
        heavy = [a for a in atom_list if a.GetAtomicNum() != 1]
        if not heavy:
            return 1.0
        nC = sum(1 for a in heavy if a.GetAtomicNum() == 6)
        return nC / len(heavy)

    # For non-ring candidates, we will remove the candidate nitrogen then see if each heavy neighbor
    # is separated into a distinct fragment. Also we use the resulting fragment as a measure of the substituent's makeup.
    for candidate in candidates:
        cand_idx = candidate.GetIdx()

        # Reject candidate if it appears to be part of an amide (attached to a carbonyl carbon)
        if is_adjacent_to_carbonyl(candidate):
            # Note: many amide nitrogens are not derived directly from NH3.
            continue

        heavy_neighbors = [nbr for nbr in candidate.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # Retrieve original indices for the heavy neighbors.
        neighbor_orig_idxs = [nbr.GetProp("orig_idx") for nbr in heavy_neighbors]

        connectivity_ok = True  # for acyclic candidate we require distinct fragments

        # For non-ring candidate, do connectivity test
        if not candidate.IsInRing():
            # Make an editable copy and remove the candidate nitrogen.
            rwmol = Chem.RWMol(mol_with_H)
            try:
                rwmol.RemoveAtom(cand_idx)
            except Exception as e:
                return False, f"Error during atom removal: {e}"
            new_mol = rwmol.GetMol()
            # Get fragments – each fragment is a tuple of atom indices (in the new_mol)
            frags = Chem.GetMolFrags(new_mol, asMols=False)
            # Build mapping: for each fragment, record the set of original indices (stored in "orig_idx")
            frag_mapping = {}
            for frag_id, frag in enumerate(frags):
                frag_mapping[frag_id] = set()
                for idx in frag:
                    atom_new = new_mol.GetAtomWithIdx(idx)
                    if atom_new.HasProp("orig_idx"):
                        frag_mapping[frag_id].add(atom_new.GetProp("orig_idx"))
            # Determine which fragment each heavy neighbor went to.
            neighbor_frag_ids = []
            for orig in neighbor_orig_idxs:
                found = False
                for frag_id, orig_set in frag_mapping.items():
                    if orig in orig_set:
                        neighbor_frag_ids.append(frag_id)
                        found = True
                        break
                if not found:
                    neighbor_frag_ids.append(-1)
            if len(set(neighbor_frag_ids)) != 3:
                connectivity_ok = False  # some substituents remain merged
        # For candidates in rings, we skip the connectivity test.
        
        if not connectivity_ok and not candidate.IsInRing():
            continue  # try next candidate

        # Now check the "hydrocarbyl" character of each substituent.
        # For non-ring candidates, we look into the entire fragment obtained after removal.
        # For ring candidates, we look only at the immediate neighbor environment.
        substituent_ok = True
        if not candidate.IsInRing():
            # Use the fragments from the connectivity test.
            for orig in neighbor_orig_idxs:
                frag_found = None
                for frag_id, orig_set in frag_mapping.items():
                    if orig in orig_set:
                        frag_found = frag_id
                        break
                if frag_found is None:
                    substituent_ok = False
                    break
                # Get fragment atoms from new_mol using the indices in that fragment.
                frag_atoms = [new_mol.GetAtomWithIdx(idx) for idx in frags[frag_found]]
                frac = carbon_fraction(frag_atoms)
                if frac < 0.5:
                    substituent_ok = False
                    break
        else:
            # For candidates in a ring, examine each direct neighbor plus its immediate heavy neighbors (excluding the candidate).
            for nbr in heavy_neighbors:
                neighbor_group = [nbr]
                for n2 in nbr.GetNeighbors():
                    if n2.GetIdx() == candidate.GetIdx() or n2.GetAtomicNum() == 1:
                        continue
                    neighbor_group.append(n2)
                frac = carbon_fraction(neighbor_group)
                if frac < 0.5:
                    substituent_ok = False
                    break

        if not substituent_ok:
            continue

        # If we reach here, candidate passed all tests.
        if candidate.IsInRing():
            return True, ("Found tertiary amine (in a ring) at atom index {} with 3 substituents that "
                          "show predominantly hydrocarbyl character."
                          .format(cand_idx))
        else:
            return True, ("Found tertiary amine at atom index {} with three substituents that become disconnected upon removal "
                          "and display predominantly hydrocarbyl character."
                          .format(cand_idx))
        
    return False, "Tertiary amine group found, but none passed all connectivity and substituent composition tests."

# Example usage (for testing):
if __name__ == "__main__":
    test_examples = [
        ("Tri-allate", "CC(C)N(C(C)C)C(=O)SCC(Cl)=C(Cl)Cl"),
        ("NK154183B", "CCC(O)CC1CCCC2(CC3OC(=O)\\C=C\\C(C)(O)C(O)C(C)C(O)C(OC4CCC(C(C)O4)N(C)C)C(O)C(C)(O)CCCCC\\C=C\\C4CC(C)(C)OC4(O)CC(O2)C3C"),
        ("triethylamine", "CCN(CC)CC"),
        ("2-[4-(dimethylamino)styryl]-1-methylpyridinium", "CN(C)c1ccc(cc1)\\C=C\\c1cccc[n+]1C"),
        ("(R)-fenpropidin", "C[C@@H](CN1CCCCC1)Cc1ccc(cc1)C(C)(C)C"),
        ("FM 1-43(2+)", "[H]C(=Cc1cc[n+](CCC[N+](CC)(CC)CC)cc1)c1ccc(cc1)N(CCCC)CCCC"),
        ("3-quinuclidinol", "OC1C[N@@]2CC[C@H]1CC2"),
        ("(R)-aceprometazine", "C[C@H](CN1c2ccccc2Sc2ccc(cc12)C(C)=O)N(C)C"),
        ("tridodecylamine", "C(CCCCCCCC)CCCN(CCCCCCCCCCCC)CCCCCCCCCCCC"),
        ("alverine", "CCN(CCCc1ccccc1)CCCc1ccccc1"),
        ("EPTC", "CCN(CCC)C(=O)SCC"),
        ("N,N-dimethylethanolamine", "CN(C)CCO")
    ]
    
    for name, smi in test_examples:
        result, reason = is_tertiary_amine(smi)
        print(f"{name}: {result} => {reason}")