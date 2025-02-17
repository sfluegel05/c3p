"""
Classifies: CHEBI:17297 UDP-sugar
"""
#!/usr/bin/env python3
"""
Classifies: UDP-sugar
Definition:
  A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to 
  an unspecified sugar via an anomeric diphosphate linkage.

This implementation uses RDKit to verify three criteria:
  1) An unmethylated uracil substructure is present.
  2) A diphosphate linkage (P-O-P) is found.
  3) There is connectivity between (at least one) phosphorus atom of a diphosphate group 
     and the uracil substructure via a short bond pathway (to ensure that the UDP moiety is intact).
"""

from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    
    A UDP-sugar must have:
      - A pyrimidine nucleotide component with an unmethylated uracil moiety.
      - A diphosphate (P-O-P) linkage.
      - A phosphorus (from the diphosphate) that is connected to the uracil substructure by a short path.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a UDP-sugar, False otherwise.
        str: Explanation of the classification.
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
    
    # Check each uracil match to ensure that none of its atoms has a terminal methyl substituent,
    # which would indicate methylation (as in thymine).
    valid_uracil_matches = []
    for match in uracil_matches:
        methyl_found = False
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider neighbors not part of the uracil ring.
                if nbr.GetIdx() not in match:
                    # If the neighbor is a carbon and is terminal (degree 1), it is a methyl.
                    if nbr.GetAtomicNum() == 6 and nbr.GetDegree() == 1:
                        methyl_found = True
                        break
            if methyl_found:
                break
        if not methyl_found:
            valid_uracil_matches.append(match)
    
    if not valid_uracil_matches:
        return False, "Only methylated uracil (thymine-like) substructures found; not UDP-sugar"
    
    # 2. Check for a diphosphate (P-O-P) linkage.
    diphosphate_smarts = "[P]-O-[P]"
    diphosphate_pattern = Chem.MolFromSmarts(diphosphate_smarts)
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    if not diphosphate_matches:
        return False, "Diphosphate (P-O-P) linkage not found"
    
    # 3. Verify connectivity between the UDP moiety and the uracil substructure.
    # Instead of checking immediate neighbors, we check for any short bond path (within a threshold)
    # between any phosphorus atom in a diphosphate motif and any atom of an unmethylated uracil.
    # This accounts for the fact that uracil is connected to a ribose which is then bonded to a phosphate.
    connectivity_found = False
    # Define a maximum bond distance threshold (e.g. 4 bonds) as acceptable for a connection.
    max_distance = 4
    
    # First, get all phosphorus atom indices that occur in a diphosphate match.
    phosphate_indices = set()
    for match in diphosphate_matches:
        # Each match returns indices corresponding to the pattern [P]-O-[P].
        # Extract phosphorus indices (atoms with atomic number 15).
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 15:
                phosphate_indices.add(idx)
    
    # Now for each phosphorus atom, check for a short linkage to any atom in any valid uracil match.
    for p_idx in phosphate_indices:
        for uracil_match in valid_uracil_matches:
            for u_idx in uracil_match:
                # Use RDKit's GetShortestPath function.
                path = Chem.GetShortestPath(mol, p_idx, u_idx)
                if path and len(path) - 1 <= max_distance:
                    connectivity_found = True
                    break
            if connectivity_found:
                break
        if connectivity_found:
            break
            
    if not connectivity_found:
        return False, "Uracil substructure not connected to a diphosphate phosphorus within a short bond distance (UDP moiety missing)"
    
    return True, "Molecule contains an unmethylated uracil substructure and a diphosphate linkage connected via a short path, consistent with a UDP-sugar"

# Example usage (for testing purposes)
if __name__ == "__main__":
    test_smiles = [
        # UDP-alpha-D-mannuronic acid
        "O[C@@H]1[C@@H](COP(O)(=O)OP(O)(=O)O[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O)C(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O",
        # UDP-4-amino-4,6-dideoxy-L-N-acetyl-beta-L-altrosamine
        "C[C@@H]1O[C@H](OP([O-])(=O)OP([O-])(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2ccc(=O)[nH]c2=O)[C@H](NC(C)=O)[C@@H](O)[C@H]1[NH3+]",
        # UDP-D-glucosamine
        "N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)OC1OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1ccc(=O)[nH]c1=O",
    ]
    
    for s in test_smiles:
        result, reason = is_UDP_sugar(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")