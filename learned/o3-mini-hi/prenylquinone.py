"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: prenylquinone

A prenylquinone is defined as a quinone substituted by a polyprenyl‐derived side chain.
That is, the molecule must contain a quinone core (i.e. a conjugated ring
that carries at least two carbonyl (C=O) groups) and
at least one substituent branch directly attached to the core that contains at least two isoprene units.
For the prenyl-derived side chain we use a SMARTS pattern for an isoprene unit "[CH2]C=C([CH3])"
and require that at least two such units are directly attached to the quinone core.
"""

from rdkit import Chem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    
    Heuristic steps:
      1. Look for a quinone core by scanning each ring in the molecule.
         A candidate quinone core is an (usually conjugated) ring (of 5 or more atoms)
         that has at least 2 carbonyl groups (C=O) attached directly to ring atoms.
      2. Look for prenyl (isoprene) units using the SMARTS pattern "[CH2]C=C([CH3])".
         Count how many of them are directly attached to one of the atoms in the prenylquinone core.
         For our purposes, we require at least two such units.
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a prenylquinone, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1: Identify a quinone core by looking at rings ---
    # We loop over all rings and check if a ring (of size >= 5) has at least 2 carbonyl substituents.
    quinone_core = None
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) < 5:
            continue  # ignore very small rings
        carbonyl_count = 0
        # For each atom in the ring, check if it has a C=O bond (i.e. double bond to an oxygen)
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider carbon atoms
            if atom.GetAtomicNum() != 6:
                continue
            # Look at each neighbor: if it is oxygen and the bond is a double bond, count it.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                        carbonyl_count += 1
                        break  # count once per ring atom
        if carbonyl_count >= 2:
            quinone_core = set(ring)
            break  # use the first ring that qualifies
    
    if quinone_core is None:
        return False, "No quinone core detected; no ring with at least two carbonyl groups attached was found"
    
    # --- Step 2: Identify prenyl (isoprene) fragments ---
    prenyl_smarts = "[CH2]C=C([CH3])"
    prenyl_pattern = Chem.MolFromSmarts(prenyl_smarts)
    if prenyl_pattern is None:
        return False, "Error creating the prenyl SMARTS pattern"
    
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern)
    if not prenyl_matches:
        return False, "No prenyl fragment detected in the molecule"
    
    # Count the number of prenyl (isoprene) units that are directly attached to the quinone core.
    attached_prenyl_count = 0
    for match in prenyl_matches:
        # For each prenyl match, check if any atom in this match is directly bonded to an atom in the quinone core.
        is_attached = False
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in quinone_core:
                    is_attached = True
                    break
            if is_attached:
                break
        if is_attached:
            attached_prenyl_count += 1

    # We require at least two attached isoprene units (i.e. a polyprenyl chain).
    if attached_prenyl_count < 2:
        return False, (f"Prenyl fragment detected but only {attached_prenyl_count} isoprene unit(s) are attached " +
                       "to the quinone core (at least 2 are required for a polyprenyl side chain)")
    
    return True, (f"Molecule contains a quinone core with {attached_prenyl_count} prenyl (isoprene) unit(s) attached")

# Example test run
if __name__ == "__main__":
    # Test with ubiquinone-2 as one example.
    test_smiles = "COC1=C(OC)C(=O)C(C\\C=C(/C)CCC=C(C)C)=C(C)C1=O"
    result, reason = is_prenylquinone(test_smiles)
    print(f"Result: {result}, Reason: {reason}")