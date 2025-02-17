"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: Tertiary amine

A tertiary amine is defined as a compound formally derived from ammonia by replacing three hydrogen atoms by hydrocarbyl groups.
This implementation first locates neutral, non‐aromatic nitrogen atoms that have exactly three non‐hydrogen neighbors.
For candidates not in a ring we then remove the candidate nitrogen and require that the three heavy substituents fall into
three distinct fragments. For each substituent (either the full fragment in acyclic cases or the immediate neighbor environment
in ring cases) we compute a “carbon fraction” (the ratio of carbon atoms to heavy atoms). In the case that the substituent is attached
through a carbonyl‐type atom (a C with at least one double bond to oxygen), we discount that atom. If all three substituents have a
carbon fraction above a threshold (here 0.4) then the tertiary amine candidate is accepted.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tertiary_amine(smiles: str):
    """
    Determines if the input molecule (given as a SMILES string) contains at least one tertiary amine group.

    A candidate tertiary amine is a neutral, non‐aromatic nitrogen that is connected to exactly three heavy atoms (i.e. non‐hydrogen).
    For non‐ring nitrogen candidates we remove the candidate and check that its three heavy neighbors fall into separate fragments.
    Then each substituent is examined to see if it displays “predominantly hydrocarbyl” character.
    (If the substituent is attached via a carbonyl carbon, that atom is discounted from the calculation.)
    
    Args:
         smiles (str): SMILES string for the molecule.
         
    Returns:
         (bool, str): A tuple where the boolean is True if a tertiary amine group is found, along with explanatory text;
                      otherwise, False and a reason why.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for correct neighbor counts.
    molH = Chem.AddHs(mol)
    # Store original atom indices in a property so we can track them after modification.
    for atom in molH.GetAtoms():
        atom.SetProp("orig_idx", str(atom.GetIdx()))
        
    # Find candidate nitrogen atoms:
    candidates = []
    for atom in molH.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        if atom.GetFormalCharge() != 0:
            continue
        if atom.GetIsAromatic():
            continue
        # We require that no explicit hydrogen is attached.
        if any(neigh.GetAtomicNum() == 1 for neigh in atom.GetNeighbors()):
            continue
        # Must have exactly three heavy (non-hydrogen) neighbors.
        heavy_nbrs = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_nbrs) != 3:
            continue
        candidates.append(atom)
    
    if not candidates:
        return False, "No candidate tertiary amine group found in the molecule."
        
    # helper: determine if an atom (assumed to be carbon) is a carbonyl carbon
    def is_carbonyl(atom):
        if atom.GetAtomicNum() != 6:
            return False
        for bond in atom.GetBonds():
            # Look for double-bonded oxygen (the bond type for a true double bond should equal 2)
            if bond.GetBondTypeAsDouble() == 2:
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8:
                    return True
        return False
    
    # helper: compute carbon fraction given a set of atoms (only heavy atoms counted).
    def carbon_fraction(atom_list):
        heavy_atoms = [a for a in atom_list if a.GetAtomicNum() != 1]
        if not heavy_atoms:
            return 1.0
        nC = sum(1 for a in heavy_atoms if a.GetAtomicNum() == 6)
        return nC / len(heavy_atoms)
    
    # For non‐ring candidates we remove the candidate N and check connectivity.
    for candidate in candidates:
        cand_idx = candidate.GetIdx()
        heavy_nbrs = [nbr for nbr in candidate.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # We also record each neighbor’s original index.
        nbr_orig_idxs = [nbr.GetProp("orig_idx") for nbr in heavy_nbrs]
        
        connectivity_ok = True
        frag_mapping = {}   # maps fragment id to set of original indices
        if not candidate.IsInRing():
            # Create an editable copy.
            rwmol = Chem.RWMol(molH)
            try:
                rwmol.RemoveAtom(cand_idx)
            except Exception as e:
                return False, f"Error during candidate removal: {e}"
            newmol = rwmol.GetMol()
            # Get fragments.
            frags = Chem.GetMolFrags(newmol, asMols=False)
            # Build mapping from each fragment (list of new indices) to the set of original indices.
            for frag in frags:
                orig_ids = set()
                for idx in frag:
                    atom_new = newmol.GetAtomWithIdx(idx)
                    if atom_new.HasProp("orig_idx"):
                        orig_ids.add(atom_new.GetProp("orig_idx"))
                # Save one fragment mapping (we simply record the set).
                frag_mapping[frozenset(frag)] = orig_ids
            # Now, for each heavy neighbor, determine which fragment (if any) contains it.
            nbr_frag_ids = []
            for orig in nbr_orig_idxs:
                found_frag = None
                for frag_set, orig_set in frag_mapping.items():
                    if orig in orig_set:
                        found_frag = frag_set
                        break
                if found_frag is None:
                    nbr_frag_ids.append(None)
                else:
                    nbr_frag_ids.append(found_frag)
            # Must have three distinct fragments.
            if None in nbr_frag_ids or len(set(nbr_frag_ids)) != 3:
                connectivity_ok = False
        
        # For substituent evaluation we will now determine an approximate carbon fraction for each group.
        substituent_ok = True
        # We relax connectivity check for ring candidates.
        # For non‐ring candidates, use the fragment obtained after removal.
        # For ring candidates, simply define the substituent to be the heavy neighbor and its immediate heavy neighbors (except candidate).
        THRESHOLD = 0.4  # minimal acceptable carbon fraction
        if not candidate.IsInRing() and connectivity_ok:
            # For each heavy neighbor of candidate, locate its fragment.
            # Build a lookup: map each neighbor (by original index) to the list of atoms in its fragment.
            for nbr, orig in zip(heavy_nbrs, nbr_orig_idxs):
                frag_atoms = None
                for frag_set, orig_set in frag_mapping.items():
                    if orig in orig_set:
                        # Collect atoms from newmol corresponding to indices in frag_set.
                        frag_atoms = [newmol.GetAtomWithIdx(i) for i in frag_set]
                        break
                if frag_atoms is None:
                    substituent_ok = False
                    break
                # Now, check if the atom attached to candidate (i.e. nbr) is a carbonyl. To do that we map original index.
                att_atom = nbr
                discount = False
                if att_atom.GetAtomicNum() == 6 and is_carbonyl(att_atom):
                    discount = True
                # If discount, remove the attachment atom from the calculation.
                if discount and len(frag_atoms) > 1:
                    # Exclude one atom (the attachment) from both numerator and denominator.
                    heavy_atoms = [a for a in frag_atoms if a.GetAtomicNum() != 1 and a.GetIdx() != att_atom.GetIdx()]
                    if heavy_atoms:
                        nC = sum(1 for a in heavy_atoms if a.GetAtomicNum() == 6)
                        frac = nC / len(heavy_atoms)
                    else:
                        frac = 1.0
                else:
                    frac = carbon_fraction(frag_atoms)
                if frac < THRESHOLD:
                    substituent_ok = False
                    break
        else:
            # For ring candidates or those failing connectivity, use a simplified evaluation:
            for nbr in heavy_nbrs:
                subgroup = [nbr]
                for n2 in nbr.GetNeighbors():
                    if n2.GetIdx() == candidate.GetIdx() or n2.GetAtomicNum() == 1:
                        continue
                    subgroup.append(n2)
                discount = False
                if nbr.GetAtomicNum() == 6 and is_carbonyl(nbr):
                    discount = True
                if discount and len(subgroup) > 1:
                    heavy_atoms = [a for a in subgroup if a.GetAtomicNum() != 1 and a.GetIdx() != nbr.GetIdx()]
                    if heavy_atoms:
                        nC = sum(1 for a in heavy_atoms if a.GetAtomicNum() == 6)
                        frac = nC / len(heavy_atoms)
                    else:
                        frac = 1.0
                else:
                    frac = carbon_fraction(subgroup)
                if frac < THRESHOLD:
                    substituent_ok = False
                    break
                    
        if not connectivity_ok and not candidate.IsInRing():
            # candidate fails connectivity test (substituents not distinct)
            continue
        if not substituent_ok:
            continue
            
        # If we reach here, candidate passes all tests.
        if candidate.IsInRing():
            return True, ("Found tertiary amine (in a ring) at atom index {} with 3 substituents that show "
                          "predominantly hydrocarbyl character."
                          .format(cand_idx))
        else:
            return True, ("Found tertiary amine at atom index {} with three substituents that, upon candidate removal, "
                          "are in distinct fragments and display predominantly hydrocarbyl character."
                          .format(cand_idx))
    
    return False, ("Tertiary amine group found, but none passed all connectivity and substituent composition tests.")


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