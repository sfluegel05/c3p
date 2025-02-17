"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
#!/usr/bin/env python3
"""
Classifies: Quaternary ammonium ion

Definition:
 A derivative of ammonium, NH4(+), in which all four hydrogens bonded to nitrogen 
 have been replaced with univalent (usually organyl) groups.

Improvements over the previous attempt:
  - Instead of extracting fragments and excluding structures with carboxylate groups,
    we directly search for candidate N+ centers.
  - We require that the candidate nitrogen is tetra‑coordinated, carries a formal +1 charge,
    has no attached hydrogens, and is non‐aromatic (typically quaternary ammonium ions are sp3).
  - For each neighbor of the candidate, we check that it is “organic” –
    namely, if it is not carbon then it must be attached (beyond the bond to N)
    to at least one carbon.
    
This approach helps retain valid quaternary ammonium ions (even when attached to other functional groups)
while reducing false positives.
"""

from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a given SMILES string contains a quaternary ammonium ion.

    A quaternary ammonium ion (QA ion) is defined as a nitrogen (N) with formal charge +1,
    exactly four bonds (tetracoordinate), and no attached (explicit or implicit) hydrogens.
    Additionally, each substituent on the nitrogen must be "organic", meaning that if the
    substituent is not a carbon atom it must be attached (besides to the candidate N) to a carbon.

    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if a qualifying quaternary ammonium ion is found, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over atoms to find candidate nitrogen atoms.
    for atom in mol.GetAtoms():
        # Check that the atom is nitrogen.
        if atom.GetAtomicNum() != 7:
            continue
        # Must have formal charge +1.
        if atom.GetFormalCharge() != 1:
            continue
        # Must have four bonds (degree = 4) and no attached hydrogens.
        if atom.GetDegree() != 4 or atom.GetTotalNumHs() != 0:
            continue
        # Optionally, require that the candidate is non‐aromatic
        if atom.GetIsAromatic():
            continue
        # (Optionally: one could also check that its hybridization is sp3)
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Now verify that each substituent (neighbor) appears to be organic.
        candidate_valid = True
        for nbr in atom.GetNeighbors():
            nbr_atomic = nbr.GetAtomicNum()
            # If neighbor is carbon, accept.
            if nbr_atomic == 6:
                continue
            # For oxygen: often oxygen substituents are allowed if they are part of an alkoxy group.
            # Here, if the neighbor is oxygen but attached (besides the candidate N) to any carbon, accept.
            if nbr_atomic == 8:
                attached_to_c = any(nbr2.GetAtomicNum() == 6 and nbr2.GetIdx() != atom.GetIdx() for nbr2 in nbr.GetNeighbors())
                if attached_to_c:
                    continue
                else:
                    candidate_valid = False
                    break
            # For any other heteroatom, require that it is connected (aside from the candidate N)
            # to at least one carbon atom.
            attached_to_c = any(nbr2.GetAtomicNum() == 6 and nbr2.GetIdx() != atom.GetIdx() for nbr2 in nbr.GetNeighbors())
            if not attached_to_c:
                candidate_valid = False
                break
        
        if candidate_valid:
            # For explanation, report the fragment size (number of heavy atoms in the whole molecule in this simplified approach).
            heavy_atom_count = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
            reason = ("Found candidate N+ with formal charge +1, 4 bonds, no hydrogens, and organic substituents. "
                      "Total heavy atoms in molecule: {}.".format(heavy_atom_count))
            return True, reason

    return False, "No qualifying quaternary ammonium ion found"


# Example usage (for testing):
if __name__ == "__main__":
    # List a few representative test cases (as provided, with notes)
    test_cases = [
        ("P(OCC[N+](C)(C)C)(OCC(OC(=O)CCCCCCCCCCCCC=1OC(CCCCC)=CC1C)COC(=O)CCCCCCCCCCC=2OC(=C(C2C)C)CCCCC)(O)=O", 
         "(3-{[11-(3,4-dimethyl-5-pentylfuran-2-yl)undecanoyl]oxy}-2-{[13-(3-methyl-5-pentylfuran-2-yl)tridecanoyl]oxy}propoxy)[2-(trimethylazaniumyl)ethoxy]phosphinic acid"),
        ("CCCCCCCCCCCCCCCC(=O)OCC(COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCC(=O)\\C=C\\C(O)=O", 
         "1-Palmitoyl-2-(5-keto-6-octendioyl)-sn-glycero-3-phosphatidylcholine"),
        ("C[N+](C)(C)CC(O)O", "Betaine aldehyde hydrate"),
        ("CC[N+](C)(CC)CCC[n+]1c(-c2ccccc2)c2cc(N)ccc2c2ccc(N)cc12", "propidium (should be detected)"),
        ("[F-][N+](F)(F)F", "Tetrafluoroammonium (should be rejected)"),
        ("C[C@H](O)[C@@H](O)COP(=O)(O)OCC[N+](C)(C)C", "A phosphocholine derivative")
    ]
    
    for smi, desc in test_cases:
        result, explanation = is_quaternary_ammonium_ion(smi)
        print("SMILES:", smi)
        print("Description:", desc)
        print("Result:", result)
        print("Explanation:", explanation)
        print("-" * 60)