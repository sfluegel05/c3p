"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: Tertiary amine

A tertiary amine is defined as a compound formally derived from ammonia by replacing
all three hydrogen atoms with hydrocarbyl groups.
This implementation first identifies all nitrogen atoms (atomic number 7)
that are neutral, non‐aromatic, and have exactly three heavy (non‐H) neighbors with no attached hydrogens.
For candidates that are not inside a ring, we then remove the candidate nitrogen (tracking each heavy neighbor via an "orig_idx")
and require that each neighbor ends up in a distinct fragment.
This connectivity test helps us avoid mis‐classifying sites that are simply embedded in a larger fused scaffold.
For candidate nitrogens in a ring we accept them directly.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tertiary_amine(smiles: str):
    """
    Determines whether the input molecule (given as a SMILES string) qualifies as a tertiary amine.
    
    Criteria (improved):
      - The molecule is parsed and explicit H's added.
      - We look for nitrogen atoms (atomic number 7) that are neutral and non‐aromatic.
      - Candidate nitrogen must have exactly three heavy (non‐H) neighbors and no attached hydrogens.
      - If the candidate is in a ring, it is accepted directly.
      - For non‐ring candidates, we remove the candidate nitrogen from a copy of the molecule and,
        tracking each heavy neighbor by an "orig_idx" property, we check that these neighbors fall into three distinct fragments.
        This ensures that the candidate’s substituents come solely off the N.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as having (at least one) tertiary amine group, False otherwise.
        str: An explanation of the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that counts are correct.
    mol_with_H = Chem.AddHs(mol)
    
    # For reliable tracking later, assign each atom a property "orig_idx" that holds its original index.
    for atom in mol_with_H.GetAtoms():
        atom.SetProp("orig_idx", str(atom.GetIdx()))
    
    candidates = []
    # Iterate over atoms looking for candidate nitrogen atoms.
    for atom in mol_with_H.GetAtoms():
        if atom.GetAtomicNum() != 7:  # Must be nitrogen.
            continue
        # Must be neutral.
        if atom.GetFormalCharge() != 0:
            continue
        # Exclude aromatic nitrogens (their lone pairs are delocalized and they are not derived directly from NH3).
        if atom.GetIsAromatic():
            continue
        # No hydrogen should be attached.
        if any(neigh.GetAtomicNum() == 1 for neigh in atom.GetNeighbors()):
            continue
        # Candidate must have exactly three heavy (non-H) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 3:
            continue
        # This atom is a candidate tertiary amine.
        candidates.append(atom)
    
    if not candidates:
        return False, "No tertiary amine group found in the molecule."
    
    # Evaluate each candidate.
    for candidate in candidates:
        cand_idx = candidate.GetIdx()
        heavy_neighbors = [nbr for nbr in candidate.GetNeighbors() if nbr.GetAtomicNum() != 1]
        # Get list of original indices for neighbors.
        neighbor_orig_idxs = [nbr.GetProp("orig_idx") for nbr in heavy_neighbors]
        
        # If candidate nitrogen is in a ring, accept it immediately.
        if candidate.IsInRing():
            return True, ("Found tertiary amine (in a ring) at atom index {} with 3 heavy substituents."
                          .format(cand_idx))
        
        # For acyclic candidate, perform connectivity test by completely removing the candidate nitrogen.
        # Create an editable copy of mol_with_H.
        rwmol = Chem.RWMol(mol_with_H)
        # Remove the candidate nitrogen.
        # Note: removal will shift atom indices. Therefore we work on a copy and use stored "orig_idx" properties.
        try:
            rwmol.RemoveAtom(cand_idx)
        except Exception as e:
            return False, "Error in atom removal: " + str(e)
        new_mol = rwmol.GetMol()
        
        # Get connected fragments (each fragment as a tuple of atom indices in new_mol).
        frags = Chem.GetMolFrags(new_mol, asMols=False)
        
        # Build a mapping from fragment id to the set of "orig_idx" values in that fragment.
        frag_mapping = {}
        for frag_id, frag in enumerate(frags):
            frag_mapping[frag_id] = set()
            for idx in frag:
                atom = new_mol.GetAtomWithIdx(idx)
                # Retrieve the original index property.
                orig = atom.GetProp("orig_idx") if atom.HasProp("orig_idx") else None
                if orig is not None:
                    frag_mapping[frag_id].add(orig)
        # Determine in which fragment each neighbor landed.
        neighbor_frag_ids = []
        for orig_idx in neighbor_orig_idxs:
            found = False
            for frag_id, orig_set in frag_mapping.items():
                if orig_idx in orig_set:
                    neighbor_frag_ids.append(frag_id)
                    found = True
                    break
            if not found:
                neighbor_frag_ids.append(-1)
        
        # If each heavy neighbor is in a different fragment,
        # then removal of the candidate disconnects the substituents.
        if len(set(neighbor_frag_ids)) == 3:
            return True, ("Found tertiary amine at atom index {} with three substituents that become disconnected upon removal."
                          .format(cand_idx))
    # If no candidate passed the connectivity test.
    return False, "Tertiary amine group found but none passed the connectivity test."

# Example usages:
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