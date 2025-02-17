"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: Tertiary amine

A tertiary amine is defined as a compound formally derived from ammonia by replacing
all three hydrogen atoms with hydrocarbyl groups. In our implementation we first look for
a neutral sp³ nitrogen that has exactly three heavy‐atom (non‐H) neighbors and no attached hydrogens.
Then, if the candidate is not in a ring we break the bonds between that N and its three neighbors and
check whether each substituent appears in a distinct disconnected fragment. This “connectivity test”
helps avoid mis‐classifying cases in which an N–substituent motif is embedded in a larger scaffold.
If the candidate is in a ring we accept it directly.
Note that some borderline cases (for example, when one substituent is “decorated” by additional groups)
might fail the connectivity test even though they are formally tertiary amines.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tertiary_amine(smiles: str):
    """
    Determines whether the input molecule (given as a SMILES string) qualifies as a tertiary amine.
    
    Criteria:
      – The molecule must contain at least one neutral (formal charge 0) sp³-hybridized nitrogen atom.
      – That nitrogen must be attached to exactly three heavy (non‐H) atoms and have no hydrogens attached.
      – If the candidate nitrogen is not in a ring, we remove the bonds connecting it to its neighbors.
         In a “true” tertiary amine each substituent (the heavy neighbor) should fall into a different disconnected fragment.
         (This test helps reject cases where the N function is embedded in a larger, fused scaffold.)
      – If the candidate nitrogen is in a ring, we accept it (since many correct tertiary amine groups lie in rings).
    
    Args:
        smiles (str): A SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a tertiary amine by our criteria, False otherwise.
        str: An explanation of the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that our counts are correct
    mol_with_H = Chem.AddHs(mol)
    
    # List to gather candidate tertiary amine nitrogen atoms
    candidates = []
    
    # Iterate over atoms looking for a candidate nitrogen.
    for atom in mol_with_H.GetAtoms():
        # Check atomic number (must be nitrogen)
        if atom.GetAtomicNum() != 7:
            continue
        # Must be neutral.
        if atom.GetFormalCharge() != 0:
            continue
        # We require that the nitrogen is sp³ hybridized.
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        # No hydrogen should be attached (explicit H, now present after AddHs)
        if any(neigh.GetAtomicNum() == 1 for neigh in atom.GetNeighbors()):
            continue
        # Candidate must have exactly three heavy (non‐H) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 3:
            continue
        # This atom is a candidate tertiary amine.
        candidates.append(atom)
    
    if not candidates:
        return False, "No tertiary amine group found in the molecule."
    
    # Now evaluate each candidate.
    for candidate in candidates:
        cand_idx = candidate.GetIdx()
        heavy_neighbors = [nbr for nbr in candidate.GetNeighbors() if nbr.GetAtomicNum() != 1]
        neighbor_indices = sorted([nbr.GetIdx() for nbr in heavy_neighbors])
        
        # If the candidate nitrogen is in a ring, we accept it as is.
        if candidate.IsInRing():
            return True, ("Found tertiary amine (in a ring) at atom index {} with 3 heavy substituents."
                          .format(cand_idx))
                          
        # For an acyclic candidate, perform a connectivity test.
        # Make a copy as an editable molecule (working on the version with explicit Hs so that bonds are clear)
        emol = Chem.EditableMol(mol_with_H)
        # Remove the bonds between candidate N and its heavy neighbors.
        for nbr_idx in neighbor_indices:
            # RemoveBond does not change atom indices.
            try:
                emol.RemoveBond(cand_idx, nbr_idx)
            except Exception:
                pass
        # Get the modified molecule after bond removals.
        mol_removed = emol.GetMol()
        # Get connected fragments (each as a tuple of atom indices)
        fragments = Chem.GetMolFrags(mol_removed, asMols=False)
        # Create a mapping from atom index to fragment id.
        frag_assignment = {}
        for frag_id, frag in enumerate(fragments):
            for idx in frag:
                frag_assignment[idx] = frag_id
        # See in which fragment each neighbor falls.
        neighbor_frags = [frag_assignment.get(idx, -1) for idx in neighbor_indices]
        # If the three heavy substituents fall into three different fragments,
        # then the candidate nitrogen serves as the sole point of connection.
        if len(set(neighbor_frags)) == 3:
            return True, ("Found tertiary amine at atom index {} with three substituents that become disconnected upon bond removal."
                          .format(cand_idx))
        # Otherwise, the candidate appears integrated in a larger scaffold.
    # If none of the candidate nitrogens tested passed our connectivity test, return False.
    return False, "Tertiary amine found but the substituents appear integrated into a larger scaffold."

# Example usage:
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