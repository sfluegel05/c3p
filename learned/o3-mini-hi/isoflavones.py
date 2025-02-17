"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: isoflavones – any isoflavonoid with a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton and its substituted derivatives.
"""

from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule belongs to the isoflavones class based on its SMILES string.
    Isoflavones are defined as isoflavonoids containing a 3-aryl-1-benzopyran-4-one scaffold.
    
    Approach:
    1. Identify a benzopyran-4-one (chromen-4-one) core using a SMARTS pattern.
    2. Search from sp2 atoms of the core for an external aromatic substituent that
       is part of a six-membered aromatic ring (as expected for the typical 3-aryl substituent).
       
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an isoflavone, False otherwise.
        str: Reason explaining the classification decision.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Define a SMARTS for a benzopyran-4-one (chromen-4-one) core.
    core_smarts = "c1ccc2c(c1)oc(=O)c(c2)"  # simple pattern to capture the fused benzene-pyran with a carbonyl.
    core = Chem.MolFromSmarts(core_smarts)
    if core is None:
        return False, "Error generating core SMARTS"
    
    core_matches = mol.GetSubstructMatches(core)
    if not core_matches:
        return False, "Molecule does not contain a benzopyran-4-one (chromen-4-one) core."
    
    # Retrieve the ring information for the whole molecule once.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Step 2. Search for a 3-aryl substituent.
    # In a typical isoflavone, one of the sp2 atoms of the chromen-4-one core will bear an external aromatic substituent.
    for match in core_matches:
        core_idx_set = set(match)
        # Iterate over atoms in the core match.
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            # Check that the atom is sp2 hybridized (typical for aromatic carbon attachment).
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP2:
                continue
            # Check neighbor atoms that are not part of the core.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in core_idx_set:
                    continue  # skip if neighbor belongs to the core
                # Consider only aromatic neighbors.
                if not nbr.GetIsAromatic():
                    continue
                # Check if the neighbor is part of at least one six-membered ring.
                for ring in ring_info:
                    if nbr.GetIdx() in ring and len(ring) == 6:
                        # Verify that every atom in this ring is aromatic.
                        if all(mol.GetAtomWithIdx(a_idx).GetIsAromatic() for a_idx in ring):
                            return True, "Molecule contains a chromen-4-one core with an external six-membered aromatic (aryl) substituent."
    
    return False, "Molecule has a benzopyran-4-one core but lacks an appropriate external aryl substituent at the 3-position."

# Example test cases; these illustrate both true and false outcomes.
if __name__ == "__main__":
    test_smiles = [
        "COc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",  # formononetin (expected True)
        "O1C2=C(CC=C(C)C)C(O)=CC(O)=C2C(=O)C(C3=CC=C(O)C=C3)=C1",  # Lupiwighteone (expected True)
        "C1=CC=C2C(=C1)C=C(C(=O)O2)C3=NN(C=C3C=C4C(=O)NC(=S)NC4=O)CCC(=O)O",  # Provided example that previously caused an error.
    ]
    for smi in test_smiles:
        result, reason = is_isoflavones(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")