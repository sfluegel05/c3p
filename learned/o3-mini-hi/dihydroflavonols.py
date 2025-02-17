"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: dihydroflavonols
Definition: Any hydroxyflavanone in which a hydroxy group is present at position 3 of the heterocyclic ring.
That is, a flavanone (benzopyran-4-one scaffold) with a saturated (dihydro) C2-C3 bond, where at least one of the saturated carbons (C3) carries an OH.
"""

from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a given SMILES string corresponds to a dihydroflavonol.
    We use an algorithmic approach: look for a 6-membered ring (the chromanone ring)
    that contains:
      - exactly one heterocyclic oxygen,
      - exactly one carbonyl carbon (C=O) (the C4 position),
      - two adjacent saturated (sp3, non‐aromatic) carbons (C2 and C3) with at least one bearing an OH group,
    and that this ring is fused to an aromatic ring (the B-ring).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a dihydroflavonol, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    found_candidate = False
    candidate_reason = ""
    
    # Utility: check if an atom has a carbonyl neighbor (a double bonded O)
    def is_carbonyl(atom):
        if atom.GetSymbol() != 'C':
            return False
        for bond in atom.GetBonds():
            # Check for a double bond to oxygen
            if bond.GetBondTypeAsDouble() == 2.0:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetSymbol() == 'O':
                    return True
        return False

    # Utility: check if a given atom (part of candidate ring) carries an -OH group out-of-ring.
    def has_exocyclic_OH(atom, ring_set):
        for nbr in atom.GetNeighbors():
            # If the neighbor is oxygen and not in the same ring, and is attached by a single bond, check if it has an implicit or explicit hydrogen.
            if nbr.GetSymbol() == 'O' and nbr.GetIdx() not in ring_set:
                # RDKit stores implicit H count.
                if nbr.GetTotalNumHs() > 0:
                    return True
        return False

    # Search for a candidate 6-membered ring having the key features
    for ring in ring_info:
        if len(ring) != 6:
            continue  # we require a 6-membered scaffold
        ring_set = set(ring)
        
        # Count heteroatoms: we require exactly one oxygen in the ring.
        oxygens_in_ring = [i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == 'O']
        if len(oxygens_in_ring) != 1:
            continue

        # Identify carbonyl carbon(s) in the ring. We require exactly one carbon which is a carbonyl (C=O).
        carbonyl_indices = [i for i in ring if is_carbonyl(mol.GetAtomWithIdx(i))]
        if len(carbonyl_indices) != 1:
            continue
        
        # Identify candidates among the remaining ring carbons for being the saturated (sp3) centers.
        # We expect two adjacent carbons that are non-aromatic and sp3 (dihydro part).
        saturated_indices = []
        for i in ring:
            atom = mol.GetAtomWithIdx(i)
            # We consider an atom "saturated" if it is not aromatic and its hybridization is sp3.
            if not atom.GetIsAromatic() and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                saturated_indices.append(i)
        if len(saturated_indices) < 2:
            continue

        # Now, we require that at least one pair of adjacent atoms in the ring (according to the ring connectivity) 
        # are both saturated and that (at least one of them) carries an exocyclic OH.
        found_saturation_adjacent = False
        for idx in range(len(ring)):
            a_idx = ring[idx]
            b_idx = ring[(idx+1) % len(ring)]  # adjacent in the ring
            if (a_idx in saturated_indices) and (b_idx in saturated_indices):
                # Check if at least one of these two carbons has an exocyclic OH
                if has_exocyclic_OH(mol.GetAtomWithIdx(a_idx), ring_set) or has_exocyclic_OH(mol.GetAtomWithIdx(b_idx), ring_set):
                    found_saturation_adjacent = True
                    break
        if not found_saturation_adjacent:
            continue

        # Finally, check for ring fusion: the flavanone scaffold should be fused to an aromatic ring (the B-ring).
        # We look for another 6-membered ring (aromatic) that shares at least 2 atoms with the candidate ring.
        fused_aromatic_found = False
        for other_ring in ring_info:
            if tuple(other_ring) == ring:
                continue
            # Check if the other ring is aromatic; a ring is aromatic if all its atoms are marked aromatic.
            if len(other_ring) == 6 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in other_ring):
                # Count shared atoms
                shared = set(ring).intersection(other_ring)
                if len(shared) >= 2:
                    fused_aromatic_found = True
                    break
        if not fused_aromatic_found:
            continue

        # If we reach here, we found a candidate dihydroflavonol scaffold.
        found_candidate = True
        candidate_reason = ("Found 6-membered chromanone-like ring with exactly one heterocyclic oxygen, one carbonyl, "
                            "adjacent saturated carbons (one with an exocyclic OH) and fused to an aromatic ring.")
        break

    if found_candidate:
        return True, candidate_reason
    else:
        return False, "Molecule does not contain the required dihydroflavonol scaffold features."

# Example usage (for testing purposes); when run as main, test with a few example SMILES.
if __name__ == "__main__":
    test_smiles = [
        "O1C(C(O)C(=O)C=2C1=CC=3OCOC3C2O)C4=CC=CC=C4",  # 3,5-Dihydroxy-6,7-methylenedioxyflavanone
        "O[C@@H]1[C@H](Oc2cc(O)ccc2C1=O)c1ccc(O)cc1",    # garbanzol
        "CC(C)=CCc1c(O)cc(O)c2C(=O)[C@H](O)[C@H](Oc12)c1ccccc1",  # glepidotin B
        "OC1C(Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1",        # (+)-taxifolin
        "COC1=C(O)C=CC(=C1)[C@H]1OC2=CC(O)=CC(O)=C2C(=O)[C@@H]1O", # dihydroisorhamnetin (example)
    ]
    for s in test_smiles:
        result, reason = is_dihydroflavonols(s)
        print(f"SMILES: {s}\nClassified as dihydroflavonol? {result}\nReason: {reason}\n")