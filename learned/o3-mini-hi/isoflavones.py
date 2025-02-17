"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: isoflavones – any isoflavonoid with a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton and its substituted derivatives.

The approach is two-step:
1. We first search for a benzopyran-4-one (chromen-4-one) fused ring system using a SMARTS pattern.
2. Then we check that one of the “free” (i.e. non-core) bonds from an sp2 carbon in that core leads to an external aromatic ring (typically the 3-aryl substituent).
If both criteria are met, we classify the molecule as an isoflavone.
"""

from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule belongs to the isoflavones class based on its SMILES string.
    Isoflavones are defined as isoflavonoids containing a 3-aryl-1-benzopyran-4-one scaffold.
    
    The algorithm first finds a benzopyran-4-one (chromen-4-one) core using a SMARTS pattern.
    Then it checks whether one of the core carbons (which is expected to have the 3-aryl substituent)
    bears an external aromatic substituent (whose ring size is 6).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an isoflavone, False otherwise.
        str: Reason explaining the classification decision.
    """
    # Parse SMILES string using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Define a SMARTS for a benzopyran-4-one (chromen-4-one) core.
    # This pattern is intended to capture a fused ring system where a benzene ring is fused to a pyran bearing a carbonyl.
    # Note: This pattern does not require the substituent yet, so it may match molecules such as coumarins.
    core_smarts = "c1ccc2c(c1)oc(=O)c(c2)"
    core = Chem.MolFromSmarts(core_smarts)
    if core is None:
        return False, "Error generating core SMARTS"
    
    core_matches = mol.GetSubstructMatches(core)
    if not core_matches:
        return False, "Molecule does not contain a benzopyran-4-one (chromen-4-one) core."
    
    # Step 2. Check for a 3-aryl substituent.
    # In the typical isoflavone, a phenyl (or substituted phenyl) ring is attached directly to the “3‐position”
    # of the chromen-4-one core. We try to identify one of the atoms in the core that has a bond to
    # an external aromatic ring (six-membered).
    # We use a benzene pattern to identify ideal candidates for the aryl substituent.
    benzene = Chem.MolFromSmarts("c1ccccc1")
    
    # For each match of the core in the molecule...
    for match in core_matches:
        core_idx_set = set(match)
        # Iterate over atoms in the core match
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            # We restrict our search to sp2 carbons (and possibly heteroatoms) in the core.
            # (For a typical isoflavone the attachment is from an sp2 carbon.)
            if atom.GetHybridization() not in (Chem.rdchem.HybridizationType.SP2,):
                continue
            # Examine neighbors of this core atom that are not part of the core.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in core_idx_set:
                    continue  # neighbor is part of core; skip
                # Check that the neighbor is aromatic and part of a six-membered ring.
                if nbr.GetIsAromatic():
                    # Count ring size(s) for the neighbor.
                    rings = nbr.GetRingInfo().NumAtomRings()
                    # (A simple check: if the neighbor is part of a benzene ring, then it will match benzene pattern.)
                    # Alternatively, we can directly check if the neighbor (or the fragment it belongs to) matches benzene.
                    # Here we use substructure matching on the neighbor's environment.
                    env = Chem.PathToSubmol(mol, Chem.rdmolops.FindAtomEnvironmentOfRadiusN(mol, 1, nbr.GetIdx()))
                    if env.HasSubstructMatch(benzene):
                        # We found an external aromatic ring attached to the core.
                        return True, "Molecule contains a chromen-4-one core with an external aromatic (aryl) substituent typical of isoflavones."

    # If no appropriate external aryl substituent was found then reject the molecule.
    return False, "Molecule has a benzopyran-4-one core but lacks an external aryl substituent at the 3-position."

# Example test cases (the user provided many examples; here are a couple)
if __name__ == "__main__":
    test_smiles = [
        "COc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",  # formononetin (expected True)
        "O1C2=C(CC=C(C)C)C(O)=CC(O)=C2C(=O)C(C3=CC=C(O)C=C3)=C1",  # Lupiwighteone (expected True)
        "COc1c(CC=C(C)C)c(O)cc2oc(=O)c(cc12)-c1ccc(O)cc1",  # Example structure that might be rejected
    ]
    for smi in test_smiles:
        result, reason = is_isoflavones(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")