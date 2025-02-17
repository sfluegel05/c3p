"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: mononitrophenol
Definition: A nitrophenol is defined as a phenol (a benzene ring having exactly one free -OH substituent) 
that carries exactly one nitro group ([N+](=O)[O-]) attached to that same ring.
"""

from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
      
    Criteria:
      1. The molecule must contain exactly one nitro group ([N+](=O)[O-]).
      2. The nitro group must be attached (through its nitrogen) to an aromatic carbon.
      3. There must be a benzene ring (6-membered aromatic ring composed solely of carbons) 
         that contains that aromatic carbon.
      4. In that benzene ring there must be exactly one "free" hydroxyl (-OH) substituent.
         By "free" we mean:
           - It is connected by a single bond to a benzene carbon.
           - The oxygen atom has atomic number 8, formal charge 0, 
             and at least one hydrogen atom (as detected by GetTotalNumHs).
           - It is not part of the nitro group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a mononitrophenol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES and add explicit hydrogens for proper H-counting.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # 1. Find exactly one nitro group.
    nitro_smarts = "[N+](=O)[O-]"
    nitro_pattern = Chem.MolFromSmarts(nitro_smarts)
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups; expected exactly one"
    
    nitro_match = nitro_matches[0]
    # Identify the nitro nitrogen (atomic number 7) in the nitro match.
    nitro_N = None
    for idx in nitro_match:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 7:
            nitro_N = atom
            break
    if nitro_N is None:
        return False, "Could not identify the nitro nitrogen atom"
    
    # 2. The nitro group must be attached to an aromatic carbon.
    attached_aromatic_c = None
    for nb in nitro_N.GetNeighbors():
        if nb.GetAtomicNum() == 6 and nb.GetIsAromatic():
            attached_aromatic_c = nb
            break
    if attached_aromatic_c is None:
        return False, "Nitro group is not attached to an aromatic carbon"
    attached_idx = attached_aromatic_c.GetIdx()
    
    # 3. Identify candidate benzene rings (6-membered, all aromatic carbons) that contain the attached aromatic carbon.
    ring_info = mol.GetRingInfo()
    candidate_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue  # not a benzene ring candidate
        atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
        # Check that every atom in the ring is a carbon and aromatic
        if not all(a.GetAtomicNum() == 6 and a.GetIsAromatic() for a in atoms_in_ring):
            continue
        if attached_idx in ring:
            candidate_rings.append(ring)
    
    if not candidate_rings:
        return False, "No benzene ring (6 aromatic carbons) containing the nitro attachment was found"
    
    # 4. For each candidate ring, count the free hydroxyl (-OH) substituents.
    # A substituent counts as a free OH if:
    #   - It is attached via a single bond to a ring carbon.
    #   - The neighboring atom is oxygen (atomic num 8) with formal charge 0.
    #   - The oxygen has at least one hydrogen.
    #   - It is not part of the nitro group.
    for ring in candidate_rings:
        oh_substituents = set()  # use a set to avoid counting the same substituent twice
        for atom_idx in ring:
            ring_atom = mol.GetAtomWithIdx(atom_idx)
            for bond in ring_atom.GetBonds():
                # Look for bonds that go outside the ring.
                nb = bond.GetOtherAtom(ring_atom)
                if nb.GetIdx() in ring:
                    continue  # neighbor is within the ring; skip
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue  # only consider single bonds for substituents
                # Exclude oxygens that are part of the nitro group.
                if nb.GetIdx() in nitro_match:
                    continue
                if nb.GetAtomicNum() != 8:
                    continue
                if nb.GetFormalCharge() != 0:
                    continue
                # Check that this oxygen has at least one hydrogen.
                if nb.GetTotalNumHs() < 1:
                    continue
                oh_substituents.add(nb.GetIdx())
        if len(oh_substituents) == 1:
            return True, "Molecule is mononitrophenol: one nitro group is attached to a benzene ring that has exactly one free hydroxyl substituent."
    
    # If no candidate ring yields exactly one free -OH, then classification fails.
    return False, "No benzene ring with exactly one free hydroxyl substituent (phenol) attached to the nitro group was found"


# Example usage for testing:
if __name__ == "__main__":
    test_molecules = [
        ("2-Methoxy-5-nitrophenol", "COc1ccc(cc1O)[N+]([O-])=O"),
        ("1-((4-hydroxy-3-nitrophenyl)acetoxy)pyrrolidine-2,5-dione", "Oc1ccc(CC(=O)ON2C(=O)CCC2=O)cc1[N+]([O-])=O"),
        ("3-nitro-L-tyrosine", "N[C@@H](Cc1ccc(O)c(c1)[N+]([O-])=O)C(O)=O"),
        ("4-nitro-m-cresol", "Cc1cc(O)ccc1[N+]([O-])=O"),
        ("2-Amino-5-nitrophenol", "Nc1ccc(cc1O)[N+]([O-])=O"),
        ("(4-hydroxy-3-iodo-5-nitrophenyl)acetic acid", "OC(=O)Cc1cc(I)c(O)c(c1)[N+]([O-])=O"),
        ("3-O-ethylentacapone", "CCOc1cc(cc(c1O)[N+]([O-])=O)\\C=C(/C#N)C(=O)N(CC)CC"),
        ("(3-bromo-4-hydroxy-5-nitrophenyl)acetic acid", "C(=O)(O)CC1=CC(=C(C(=C1)Br)O)[N+]([O-])=O"),
        ("roxarsone (III)", "C=1C([As](O)O)=CC([N+](=O)[O-])=C(C1)[O-]"),
        ("2-hydroxy-5-nitrophenyl hydrogen sulfate", "C1(=CC=C(C=C1OS(O)(=O)=O)[N+]([O-])=O)O"),
        ("6-[(4-hydroxy-3-nitrophenyl)acetamido]caproic acid", "OC(=O)CCCCCNC(=O)Cc1ccc(O)c(c1)[N+]([O-])=O"),
        ("3-nitrotyramine", "NCCc1ccc(O)c(c1)[N+]([O-])=O"),
        ("6-[(3-bromo-4-hydroxy-5-nitrophenyl)acetamido]caproic acid", "OC(=O)CCCCCNC(=O)Cc1cc(Br)c(O)c(c1)[N+]([O-])=O"),
        ("3-nitrotyrosine", "NC(Cc1ccc(O)c(c1)[N+]([O-])=O)C(O)=O"),
        ("3-nitrophenol", "Oc1cccc(c1)[N+]([O-])=O"),
        ("6-(4-hydroxy-5-iodo-3-nitrobenzamido)hexanoic acid", "OC(=O)CCCCCNC(=O)c1cc(I)c(O)c(c1)[N+]([O-])=O"),
        ("1-((4-hydroxy-5-iodo-3-nitrophenyl)acetoxy)pyrrolidine-2,5-dione", "Oc1c(I)cc(CC(=O)ON2C(=O)CCC2=O)cc1[N+]([O-])=O"),
        ("5-nitrovanillin", "[H]C(=O)c1cc(OC)c(O)c(c1)[N+]([O-])=O"),
        ("roxarsone", "Oc1ccc(cc1[N+]([O-])=O)[As](O)(O)=O"),
        ("entacapone", "CCN(CC)C(=O)C(\\C#N)=C\\c1cc(O)c(O)c(c1)[N+]([O-])=O"),
        ("3-iodo-5-nitrophenol", "Oc1cc(I)cc(c1)[N+]([O-])=O"),
        ("2-nitrophenol", "Oc1ccccc1[N+]([O-])=O"),
        ("(4-hydroxy-3-nitrophenyl)acetyl azide", "Oc1ccc(CC(=O)N=[N+]=[N-])cc1[N+]([O-])=O"),
        ("4-nitrophenol", "Oc1ccc(cc1)[N+]([O-])=O"),
        ("tolcapone", "Cc1ccc(cc1)C(=O)c1cc(O)c(O)c(c1)[N+]([O-])=O"),
        ("2-Amino-4-nitrophenol", "Nc1cc(ccc1O)[N+]([O-])=O"),
        ("6-(4-hydroxy-3-nitrobenzamido)hexanoic acid", "OC(=O)CCCCCNC(=O)c1ccc(O)c(c1)[N+]([O-])=O"),
        ("(3-bromo-4-hydroxy-5-nitrophenyl)acetyl azide", "Oc1c(Br)cc(CC(=O)N=[N+]=[N-])cc1[N+]([O-])=O"),
        ("2-Methoxy-4-nitrophenol", "COc1cc(ccc1O)[N+]([O-])=O"),
        ("3-O-methylentacapone", "CCN(CC)C(=O)C(\\C#N)=C\\c1cc(OC)c(O)c(c1)[N+]([O-])=O")
    ]
    
    for name, smi in test_molecules:
        result, reason = is_mononitrophenol(smi)
        print(f"Test: {name}  SMILES: {smi}")
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 70)