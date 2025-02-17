"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: mononitrophenol
Definition: A mononitrophenol is defined as a phenol (a benzene ring having exactly one neutral -OH substituent)
that carries exactly one nitro group ([N+](=O)[O-]) attached to the same ring.
"""

from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    
    Criteria:
      1. The molecule must contain exactly one nitro group ([N+](=O)[O-]). 
         Additionally, that nitro group must be attached (through its N) to an aromatic carbon.
      2. There must be an aromatic benzene ring (6-membered ring, all carbons, aromatic)
         that contains the ring atom to which the nitro group is attached.
      3. That benzene ring must have exactly one “free” hydroxyl (-OH) substituent. 
         A substituent counts as a hydroxyl if:
            - It is attached by a single bond to a ring carbon.
            - It is an oxygen atom with formal charge 0.
            - It has at least one hydrogen (using GetTotalNumHs).
         (Any oxygen that is part of the nitro group is not counted.)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # 1. Identify nitro groups and require exactly one.
    nitro_smarts = "[N+](=O)[O-]"
    nitro_pattern = Chem.MolFromSmarts(nitro_smarts)
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups; expected exactly one"
    
    nitro_match = nitro_matches[0]
    # Identify the nitro nitrogen atom (atomic num 7)
    nitro_N = None
    for idx in nitro_match:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 7:
            nitro_N = atom
            break
    if nitro_N is None:
        return False, "Could not identify the nitro nitrogen atom"
    
    # 2. Identify an aromatic carbon attached to the nitro group.
    attached_aromatic_c = None
    for nb in nitro_N.GetNeighbors():
        # The nitro group should be attached to an aromatic carbon.
        if nb.GetAtomicNum() == 6 and nb.GetIsAromatic():
            attached_aromatic_c = nb
            break
    if attached_aromatic_c is None:
        return False, "Nitro group is not attached to an aromatic carbon"
    attached_idx = attached_aromatic_c.GetIdx()
    
    # 3. Get candidate benzene rings from the molecule.
    # We restrict to 6-membered rings that are aromatic and contain only carbons.
    ring_info = mol.GetRingInfo()
    candidate_rings = []
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue  # not benzene
        atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
        # Check that all atoms are carbon and aromatic.
        if not all(a.GetAtomicNum() == 6 and a.GetIsAromatic() for a in atoms_in_ring):
            continue
        if attached_idx in ring:
            candidate_rings.append(ring)
            
    if not candidate_rings:
        return False, "No candidate benzene ring (6 aromatic carbons) containing the nitro attachment was found"
    
    # 4. For each candidate ring, count free hydroxyl (-OH) substituents.
    # A substituent is considered as an -OH if:
    #  - It is attached via a single bond to a ring atom (i.e., the bond must be SINGLE).
    #  - The attached atom is oxygen (atomic number 8) with formal charge 0.
    #  - The oxygen has at least one hydrogen (using GetTotalNumHs).
    #  - It is not part of the nitro group (its idx should not be in nitro_match).
    for ring in candidate_rings:
        oh_count = 0  # Count of free –OH substituents on the ring.
        for atom_idx in ring:
            ring_atom = mol.GetAtomWithIdx(atom_idx)
            for bond in ring_atom.GetBonds():
                # We need to consider bonds from the ring atom to an external substituent.
                nb = bond.GetOtherAtom(ring_atom)
                if nb.GetIdx() in ring:
                    continue  # still inside the ring
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    continue  # only single bonds count for substituents
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
                oh_count += 1
        if oh_count == 1:
            return True, "Molecule contains one nitro group attached to a benzene ring with exactly one free hydroxyl substituent"
    
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