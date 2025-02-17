"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
#!/usr/bin/env python3
"""
Classifies: Quaternary ammonium ion

Definition:
 A derivative of ammonium, NH4(+), in which all four hydrogens bonded to nitrogen
 have been replaced with univalent (usually organyl) groups.

Approach:
 1. Parse the molecule.
 2. Look for nitrogen atoms (atomic number 7) with formal charge +1, exactly 4 bonds,
    and zero attached (explicit or implicit) hydrogens.
 3. For each candidate, check each substituent:
      - If the neighbor is carbon, accept.
      - If it is oxygen, allow it if there is at least one hydrogen (i.e. a hydroxyl group).
      - Otherwise, require that the substituent is “organic” – meaning it is attached (besides the candidate N)
        to at least one carbon.
 4. Identify the connected fragment (using GetMolFrags) that contains the candidate.
      - To avoid errors from PathToSubmol, we sort the fragment’s indices and use a try/except block.
      - If the fragment does not contain phosphorus (atomic number 15)
        and does contain a carboxylate group ([CX3](=O)[O-]), skip the candidate
        (it is likely part of a carnitine derivative or similar).
 5. If a candidate passes all checks, return True with some details.

If no candidate qualifies, return False and a reason.
"""

from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determine whether the input SMILES string represents a molecule
    that contains a quaternary ammonium ion (N+ with four organic substituents).

    Args:
      smiles (str): SMILES string for the molecule.
      
    Returns:
      bool: True if a qualifying quaternary ammonium ion is found, False otherwise.
      str: Explanation of the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all fragments (each as a tuple of atom indices).
    frags = Chem.GetMolFrags(mol, asMols=False)
    
    # SMARTS pattern for detecting a carboxylate group.
    carboxylate_smarts = Chem.MolFromSmarts("[CX3](=O)[O-]")
    
    # Loop over atoms looking for candidate nitrogen atoms.
    for atom in mol.GetAtoms():
        # Candidate must be nitrogen with formal charge +1.
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1:
            continue
        # For a proper quaternary ammonium ion, N must have exactly 4 bonds
        # and no attached explicit/implicit hydrogens.
        if atom.GetDegree() != 4 or atom.GetTotalNumHs() != 0:
            continue
        
        valid_substituents = True
        # Check each neighbor (substituent) of the candidate N.
        for nbr in atom.GetNeighbors():
            # Case 1: Accept if neighbor is carbon.
            if nbr.GetAtomicNum() == 6:
                continue
            # Case 2: Allow oxygen if it carries at least one hydrogen.
            if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() >= 1:
                continue
            # Otherwise, if the neighbor is not carbon or allowed oxygen, require that it is attached
            # (besides the candidate N) to at least one carbon.
            has_carbon = False
            for sub_nbr in nbr.GetNeighbors():
                if sub_nbr.GetIdx() == atom.GetIdx():
                    continue
                if sub_nbr.GetAtomicNum() == 6:
                    has_carbon = True
                    break
            if not has_carbon:
                valid_substituents = False
                break
        
        if not valid_substituents:
            continue
        
        # Identify which fragment (connected component) contains the candidate atom.
        frag_indices = None
        for frag in frags:
            if atom.GetIdx() in frag:
                frag_indices = frag
                break
        if frag_indices is None:
            continue  # Should not occur
        
        # Try to create a submolecule from the fragment.
        try:
            # Sort the fragment indices to help avoid indexing errors.
            sorted_indices = sorted(list(frag_indices))
            frag_mol = Chem.PathToSubmol(mol, sorted_indices)
        except Exception as e:
            # If there is an error creating the fragment molecule, skip this candidate.
            continue
        
        # Check for phosphorus (atomic number 15) in the fragment.
        has_phosphorus = any(a.GetAtomicNum() == 15 for a in frag_mol.GetAtoms())
        
        # If the fragment does NOT contain phosphorus and has a carboxylate,
        # then it is likely part of a carnitine-like structure; skip the candidate.
        if not has_phosphorus and frag_mol.HasSubstructMatch(carboxylate_smarts):
            continue
        
        # Candidate N+ qualifies.
        frag_size = frag_mol.GetNumHeavyAtoms()
        reason = ("Found candidate N+ with formal charge +1, exactly 4 bonds, and no hydrogens; "
                  "all substituents are organic. Fragment size: {} heavy atoms {}."
                  .format(frag_size, "(contains phosphorus)" if has_phosphorus else ""))
        return True, reason
    
    # No candidate was found; return a reason.
    return False, "No qualifying quaternary ammonium ion found"


# Example usage (for testing):
if __name__ == "__main__":
    test_cases = [
        ("P(OCC[N+](C)(C)C)(OCC(OC(=O)CCCCCCCCCCCCC=1OC(CCCCC)=CC1C)COC(=O)CCCCCCCCCCC=2OC(=C(C2C)C)CCCCC)(O)=O",
         "Phosphinic acid with a choline-like group"),
        ("CCCCCCCCCCCCCCCC(=O)OCC(COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCC(=O)\\C=C\\C(O)=O",
         "1-Palmitoyl-2-(5-keto-6-octendioyl)-sn-glycero-3-phosphatidylcholine"),
        ("C[N+](C)(C)CC(O)O",
         "Betaine aldehyde hydrate"),
        ("C[N+](C)(C)O",
         "Hydroxytrimethylaminium (previous false negative)"),
        # A known carnitine derivative should be skipped.
        ("CCCC\\C=C\\C\\C=C\\CCCCCCCC(=O)OC(CC([O-])=O)C[N+](C)(C)C",
         "Carnitine derivative (should not be classified as quaternary ammonium ion)"),
        # Tetrafluoroammonium (should be rejected)
        ("[F-][N+](F)(F)F",
         "Tetrafluoroammonium (should be rejected)")
    ]
    
    for smi, desc in test_cases:
        result, explanation = is_quaternary_ammonium_ion(smi)
        print("SMILES:", smi)
        print("Description:", desc)
        print("Result:", result)
        print("Explanation:", explanation)
        print("-" * 60)