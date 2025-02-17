"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
#!/usr/bin/env python3
"""
Classifies: Quaternary ammonium ion

Definition:
  A derivative of ammonium, NH4(+), in which all four hydrogens bonded to nitrogen
  have been replaced with univalent (usually organyl) groups.
  
Improved approach:
  1. Parse the molecule.
  2. Search for any nitrogen (atomic number 7) that has formal charge +1, exactly 4 bonds, and no hydrogens.
  3. For each substituent of candidate N:
       - Allow if it is carbon (atomic number 6).
       - Allow if it is oxygen (atomic number 8) and has at least 1 total hydrogen (accepting a hydroxyl group).
       - Otherwise, require that the substituent is “organic”, that is, it is attached (besides the candidate N)
         to at least one carbon.
  4. Identify the connected fragment (set of atoms) in which the candidate lies. If the fragment does NOT contain a phosphorus atom
     and it contains a carboxylate group (SMARTS “[CX3](=O)[O-]”), then we assume that it comes from a carnitine derivative 
     (or similar) and skip the candidate.
     
If a candidate obeys all these checks, we return True along with details.
If none are found, we return False and a reason.
"""

from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determine whether the input SMILES string represents a molecule
    that contains a quaternary ammonium ion (N+ with four organic substituents). 
    Organic substituents are defined as either carbons, or – for non-carbons – groups that have at least one carbon,
    with the special case of allowing hydroxyl groups (i.e. oxygen with at least one H).
    
    Additionally, if the candidate N+ is in a small fragment that also contains a carboxylate group 
    (and does not contain phosphorus), then it is considered likely to be part of a carnitine-like structure 
    and is not accepted.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if a qualifying quaternary ammonium ion was found, False otherwise.
        str: A reason explaining the decision.
              If an error occurs the function may return (None, None).
    """
    # Parse the SMILES and return error if invalid.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all fragments (each as a tuple of atom indices)
    frags = Chem.GetMolFrags(mol, asMols=False)
    
    # SMARTS pattern for detecting a carboxylate group.
    carboxylate_smarts = Chem.MolFromSmarts("[CX3](=O)[O-]")
    
    # Iterate over atoms and select candidate nitrogen atoms.
    for atom in mol.GetAtoms():
        # Must be nitrogen with formal charge +1.
        if atom.GetAtomicNum() != 7 or atom.GetFormalCharge() != 1:
            continue
        # For a proper quaternary ammonium ion, the nitrogen should have exactly 4 bonds
        # and no attached explicit/implicit hydrogens.
        if atom.GetDegree() != 4 or atom.GetTotalNumHs() != 0:
            continue
        
        valid_substituents = True
        # Check each neighbor (substituent).
        for nbr in atom.GetNeighbors():
            # Case 1: If the neighbor is a carbon, accept.
            if nbr.GetAtomicNum() == 6:
                continue
            # Case 2: Allow oxygen if it carries at least one hydrogen (implying a hydroxyl group).
            if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() >= 1:
                continue
            # Otherwise, for non-carbon/neither allowed oxygen, require that the neighbor is attached
            # (besides the candidate N) to at least one carbon.
            # (We iterate over the neighbors of the substituent.)
            has_carbon = False
            for sub_nbr in nbr.GetNeighbors():
                # Skip the candidate nitrogen itself.
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
        
        # Determine which fragment (connected component) contains this candidate atom.
        frag_indices = None
        for frag in frags:
            if atom.GetIdx() in frag:
                frag_indices = frag
                break
        if frag_indices is None:
            continue  # Should not happen
        
        # Create the sub-molecule corresponding to the fragment.
        frag_mol = Chem.PathToSubmol(mol, list(frag_indices))
        
        # Check if the fragment contains any phosphorus.
        has_P = any(a.GetAtomicNum() == 15 for a in frag_mol.GetAtoms())
        
        # If the fragment does NOT contain phosphorus and
        # if it contains a carboxylate group, then reject this candidate.
        if not has_P and frag_mol.HasSubstructMatch(carboxylate_smarts):
            # Likely part of a carnitine derivative or similar; skip.
            continue
        
        # If we reached here, candidate N+ qualifies.
        frag_size = frag_mol.GetNumHeavyAtoms()
        reason = ("Found a candidate nitrogen atom (N+) with formal charge +1, exactly 4 bonds and no attached hydrogens, "
                  "with all substituents organic. Fragment size = {} heavy atoms{}."
                  .format(frag_size, " (contains phosphorus)" if has_P else ""))
        return True, reason
    
    # No candidate was found.
    return False, "No quaternary ammonium group (N+ with 4 organic substituents) found"


# Example usage (for testing):
if __name__ == "__main__":
    test_cases = [
        # True positives:
        ("P(OCC[N+](C)(C)C)(OCC(OC(=O)CCCCCCCCCCCCC=1OC(CCCCC)=CC1C)COC(=O)CCCCCCCCCCC=2OC(=C(C2C)C)CCCCC)(O)=O",
         "Phosphinic acid with a choline-like group"),
        ("CCCCCCCCCCCCCCCC(=O)OCC(COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCC(=O)\\C=C\\C(O)=O",
         "1-Palmitoyl-2-(5-keto-6-octendioyl)-sn-glycero-3-phosphatidylcholine"),
        ("C[N+](C)(C)CC(O)O",
         "Betaine aldehyde hydrate"),
        ("C[N+](C)(C)O",
         "Hydroxytrimethylaminium (previous false negative)"),
        # A known false positive (carnitine derivative) should now be filtered out.
        ("CCCC\\C=C\\C\\C=C\\CCCCCCCC(=O)OC(CC([O-])=O)C[N+](C)(C)C",
         "Carnitine derivative (should not be classified as quaternary ammonium ion)"),
        # Additional test: tetrafluoroammonium should be rejected because substituents are non‐organic.
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