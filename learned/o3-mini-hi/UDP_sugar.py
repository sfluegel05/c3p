"""
Classifies: CHEBI:17297 UDP-sugar
"""
#!/usr/bin/env python3
"""
Classifies: UDP-sugar
Definition:
  A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to 
  an unspecified sugar via an anomeric diphosphate linkage.
  
This implementation uses RDKit to check:
  1) That an unmethylated uracil moiety is present. We first match the uracil ring identified 
     by the SMARTS "n1ccc(=O)[nH]c1=O" and then for each match verify that none of the ring atoms 
     have an external methyl substituent (a carbon neighbor having degree 1).
  2) That a diphosphate connectivity – a P-O-P bond – exists.
  3) That at least one phosphorus atom from a diphosphate bond is connected (via any bond) 
     to the uracil substructure (to "anchor" a UDP moiety), so that extra phosphates (e.g. on the sugar)
     will not cause a false negative.
If these conditions are met, the molecule is classified as a UDP-sugar.
"""

from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    
    A UDP-sugar must have:
      - A pyrimidine nucleotide component with a uracil moiety that is not methylated (i.e. not thymine).
      - A diphosphate linkage (P-O-P), meaning that at least one pair of phosphorus atoms is directly connected via an oxygen.
      - Connectivity between the uracil and one of the phosphorus atoms from the diphosphate (to verify a UDP moiety).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a UDP-sugar, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Look for an unmethylated uracil substructure.
    # SMARTS for a uracil moiety.
    uracil_smarts = "n1ccc(=O)[nH]c1=O"
    uracil_pattern = Chem.MolFromSmarts(uracil_smarts)
    uracil_matches = mol.GetSubstructMatches(uracil_pattern)
    if not uracil_matches:
        return False, "Uracil substructure not found"
    
    # Check each matched uracil ring to ensure no out-of-ring atom is a terminal methyl group.
    # (We assume that if any ring atom has an external carbon neighbor that is terminal, the ring is methylated.)
    valid_uracil_found = False
    for match in uracil_matches:
        methyl_found = False
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider atoms not in the uracil ring.
                if nbr.GetIdx() not in match:
                    # If the neighbor is carbon, check if it is a terminal methyl (degree 1).
                    if nbr.GetAtomicNum() == 6 and nbr.GetDegree() == 1:
                        methyl_found = True
                        break
            if methyl_found:
                break
        if not methyl_found:
            valid_uracil_found = True
            break
    if not valid_uracil_found:
        return False, "Methylated uracil (thymine-like) substructure found; not UDP-sugar"
    
    # 2. Check for a diphosphate (P-O-P) linkage.
    diphosphate_smarts = "[P]-O-[P]"
    diphosphate_pattern = Chem.MolFromSmarts(diphosphate_smarts)
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "Diphosphate (P-O-P) linkage not found"
    
    # 3. Verify that at least one phosphorus of a diphosphate bond is connected (by any bond) to the uracil substructure.
    linked_to_uracil = False
    for match in uracil_matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 15:  # phosphorus
                    linked_to_uracil = True
                    break
            if linked_to_uracil:
                break
        if linked_to_uracil:
            break
    if not linked_to_uracil:
        return False, "Uracil substructure not connected to any phosphorus atom (UDP moiety missing)"
    
    # Passed all checks.
    return True, "Molecule contains an unmethylated uracil substructure and a diphosphate linkage indicative of a UDP-sugar"

# Example usage (for testing purposes)
if __name__ == "__main__":
    test_smiles = [
        # UDP-alpha-D-mannuronic acid (expected True)
        "O[C@@H]1[C@@H](COP(O)(=O)OP(O)(=O)O[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O)C(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O",
        # dTDP-4-acetamido-4,6-dideoxy-alpha-D-galactose (expected False due to methyl substitution)
        "C[C@H]1O[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H](C[C@@H]2O)n2cc(C)c(=O)[nH]c2=O)[C@H](O)[C@@H](O)[C@H]1NC(C)=O",
        # UDP-N-acetyl-alpha-D-glucosamine 3-phosphate (expected True, as extra phosphate on the sugar is allowed)
        "CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](COP(O)(=O)O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O)[C@H]2O)OC1OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1ccc(=O)[nH]c1=O",
    ]
    
    for s in test_smiles:
        result, reason = is_UDP_sugar(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")