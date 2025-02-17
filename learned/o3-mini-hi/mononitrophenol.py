"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: mononitrophenol
Definition: A mononitrophenol is defined as a phenol (a benzene ring having exactly one neutral -OH substituent)
that carries exactly one nitro group ([N+](=O)[O-]) attached on the ring.
Examples include compounds such as 2-Methoxy-5-nitrophenol, 3-nitrotyrosine etc.
"""

from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.

    A mononitrophenol must:
      1. Contain exactly one nitro group ([N+](=O)[O-]).
      2. Have that nitro group attached to an aromatic carbon.
      3. Belong to a benzene ring (6-membered, all atoms being carbon) that carries exactly one neutral hydroxyl (-OH) substituent.
         (A substituent is considered a hydroxyl if the oxygen (atomic number 8) is not part of the nitro group,
          has at least one hydrogen (using GetTotalNumHs) and a formal charge of 0.)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that the hydroxyl groups have explicit hydrogen atoms.
    mol = Chem.AddHs(mol)
    
    # 1. Identify nitro groups using SMARTS. Require exactly one.
    nitro_smarts = "[N+](=O)[O-]"
    nitro_pattern = Chem.MolFromSmarts(nitro_smarts)
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups; expected exactly one"
    
    # 2. Find the aromatic carbon that is attached to the nitro group.
    nitro_match = nitro_matches[0]
    nitro_N = None  # the nitrogen atom in the nitro group
    for idx in nitro_match:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 7:  # nitrogen
            nitro_N = atom
            break
    if nitro_N is None:
        return False, "Could not identify the nitro nitrogen atom"

    # Find an aromatic carbon neighbor of the nitro nitrogen.
    attached_aromatic_c = None
    for nb in nitro_N.GetNeighbors():
        if nb.GetAtomicNum() == 6 and nb.GetIsAromatic():
            attached_aromatic_c = nb
            break
    if attached_aromatic_c is None:
        return False, "Nitro group is not attached to an aromatic carbon"
    attached_idx = attached_aromatic_c.GetIdx()
    
    # 3. Get candidate benzene rings that include the attached aromatic carbon.
    ring_info = mol.GetRingInfo()
    candidate_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) != 6:  # benzene rings are 6-membered
            continue
        atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
        # A benzene ring: all atoms should be carbon.
        if not all(a.GetAtomicNum() == 6 for a in atoms_in_ring):
            continue
        # Check that the ring contains the aromatic carbon from the nitro group.
        if attached_idx in ring:
            candidate_rings.append(ring)
    
    if not candidate_rings:
        return False, "No candidate benzene ring (6-membered, all carbon) containing the nitro attachment was found"
    
    # 4. For each candidate ring, count the hydroxyl (-OH) substituents.
    # We count a substituent if:
    #   - It is attached to a ring atom (i.e., not in the ring).
    #   - It is an oxygen atom (atomic number 8) not part of the nitro group.
    #   - It has at least one hydrogen (GetTotalNumHs() >= 1)
    #   - It has formal charge 0.
    for ring in candidate_rings:
        oh_count = 0
        # Iterate over atoms in the ring.
        for atom_idx in ring:
            ring_atom = mol.GetAtomWithIdx(atom_idx)
            # Look at each neighbor that is not in the ring.
            for nb in ring_atom.GetNeighbors():
                if nb.GetIdx() in ring:
                    continue  # neighbor is inside the ring
                # Skip if the neighbor atom is part of the nitro group.
                if nb.GetIdx() in nitro_match:
                    continue
                # Check if neighbor is oxygen.
                if nb.GetAtomicNum() != 8:
                    continue
                # Check that this oxygen is neutral (formal charge 0) and has at least one explicit hydrogen.
                if nb.GetFormalCharge() != 0:
                    continue
                # Sometimes the hydrogen may not be visible, so we use GetTotalNumHs().
                if nb.GetTotalNumHs() < 1:
                    continue
                oh_count += 1
        
        # Require exactly one phenolic -OH substituent on that ring.
        if oh_count == 1:
            return True, "Molecule contains one nitro group attached to a benzene ring with exactly one neutral hydroxyl substituent"
    
    return False, "No benzene ring with exactly one neutral hydroxyl substituent (phenol) attached to the nitro group was found"


# Example usage for testing:
if __name__ == "__main__":
    # Test with several SMILES strings:
    test_molecules = [
        ("2-Methoxy-5-nitrophenol", "COc1ccc(cc1O)[N+]([O-])=O"),
        ("3-nitro-L-tyrosine", "N[C@@H](Cc1ccc(O)c(c1)[N+]([O-])=O)C(O)=O"),
        ("4-nitro-m-cresol", "Cc1cc(O)ccc1[N+]([O-])=O"),
        ("2-Amino-5-nitrophenol", "Nc1ccc(cc1O)[N+]([O-])=O"),
        ("3-nitrophenol", "Oc1cccc(c1)[N+]([O-])=O"),
        ("4-nitrophenol", "Oc1ccc(cc1)[N+]([O-])=O")
    ]
    
    for name, smiles in test_molecules:
        result, reason = is_mononitrophenol(smiles)
        print(f"Test: {name}  SMILES: {smiles}")
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 70)