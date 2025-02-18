"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: 3-oxo-fatty acyl-CoA
An oxo fatty acyl-CoA results from the condensation of the thiol group of Coenzyme A
with the carboxy group of any 3-oxo fatty acid.
This program first checks for a Coenzyme A adenine fragment, then searches for a 3-oxo fatty acyl thioester fragment.
Two alternate SMARTS patterns are used to capture the two orientations:
   Pattern 1: R–C(=O)–CH2–C(=O)S
   Pattern 2: S–C(=O)–CH2–C(=O)–R
Additionally, to avoid false positives we verify that the methylene (CH2) atom has no extra substituents
and that the sulfur seen in the match only has the two expected connections.
"""

from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    
    Criteria:
      1. The SMILES must be valid and should not show deprotonated oxygens (e.g. "[O-]").
      2. The molecule must contain a Coenzyme A adenine fragment.
      3. The molecule must contain an uncompromised 3-oxo fatty acyl thioester fragment.
         We require a substructure matching one of two SMARTS patterns:
             a) "[#6]-C(=O)[CH2]C(=O)S" (e.g. acetoacetyl-CoA type)
             b) "S-C(=O)[CH2]C(=O)-[#6]" (reverse orientation)
         and then we further check that the bridging CH2 (the [CH2] in the pattern)
         is not substituted (i.e. it only connects to its two neighbors in the fragment)
         and that the S atom in the pattern has exactly two neighbors.
    
    Args:
      smiles (str): SMILES string representing the molecule.
      
    Returns:
      bool: True if the molecule meets the criteria, False otherwise.
      str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # We require neutral compounds (no deprotonated oxygens)
    if "[O-]" in smiles:
        return False, "Contains deprotonated oxygens; expected a neutral CoA derivative"
    
    # Check for Coenzyme A adenine fragment (common substructure in CoA)
    coa_pattern = Chem.MolFromSmarts("n1cnc2ncnc(c12)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Does not contain a Coenzyme A adenine fragment"
    
    # Define SMARTS patterns for the 3-oxo fatty acyl thioester fragment.
    # Pattern 1: R–C(=O)–CH2–C(=O)S
    frag_smarts1 = "[#6]-C(=O)[CH2]C(=O)S"
    # Pattern 2: S–C(=O)–CH2–C(=O)–R
    frag_smarts2 = "S-C(=O)[CH2]C(=O)-[#6]"
    frag_patterns = [frag_smarts1, frag_smarts2]
    
    valid_fragment_found = False
    
    for smarts in frag_patterns:
        frag_pattern = Chem.MolFromSmarts(smarts)
        if frag_pattern is None:
            continue  # Safety check
        matches = mol.GetSubstructMatches(frag_pattern)
        if not matches:
            continue
        # Check each match for additional connectivity restrictions:
        for match in matches:
            # Expecting match to return 5 atoms:
            # For pattern 1: [0]: R-group carbon, [1]: first carbonyl C, [2]: bridging CH2, [3]: second carbonyl C, [4]: thioester S.
            # For pattern 2 the order will be reversed for the S and R side.
            # Locate the CH2 atom in the match: it is the one that matches [CH2] in each pattern.
            # Searching for an atom with atomic number 6 and exactly 2 attached heavy atoms.
            ch2_ok = False
            s_ok = False
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6 and atom.GetSymbol() == "C" and atom.GetTotalNumHs() == 2:
                    # Check that its connections (neighbors excluding hydrogens) equal 2.
                    heavy_neighbors = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() != 1]
                    if len(heavy_neighbors) == 2:
                        ch2_ok = True
                        ch2_idx = idx
                        break
            # Similarly, find the sulfur atom in the match (atomic number 16) and check that it has exactly 2 neighbors.
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 16:
                    # We require that S has exactly 2 neighbors:
                    heavy_neighbors = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() != 1]
                    if len(heavy_neighbors) == 2:
                        s_ok = True
                        break
            if ch2_ok and s_ok:
                valid_fragment_found = True
                break  # No need to test further matches for this pattern.
        if valid_fragment_found:
            break
    
    if not valid_fragment_found:
        return False, "Does not contain an uncompromised 3-oxo fatty acyl thioester substructure"
    
    return True, "Contains a valid 3-oxo fatty acyl thioester fragment linked to a Coenzyme A moiety"

# Example usage:
if __name__ == "__main__":
    # Use one of the provided sample SMILES: 3-oxoadipyl-CoA.
    example_smiles = ("CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)"
                      "[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CC(=O)CCC(O)=O")
    result, reason = is_3_oxo_fatty_acyl_CoA(example_smiles)
    print(result, reason)