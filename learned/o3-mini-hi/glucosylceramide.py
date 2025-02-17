"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: glucosylceramide
Definition: any cerebroside in which the monosaccharide head group is glucose.
A true glucosylceramide should contain exactly one glucose ring (typically attached via an oxygen),
an amide bond connecting a fatty acyl chain to a sphingoid base, and a high enough molecular weight.
Furthermore, the amide bond should be “ceramide‐like”: the acyl group and the sphingoid part should
each exhibit a contiguous carbon chain of at least a specified length.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    
    Checks:
      1) The molecule contains exactly one glucose ring.
      2) The molecule contains at least one amide bond.
      3) At least one amide bond is "ceramide‐like": the carbonyl carbon is attached to a long (fatty acyl) chain 
         and the amide nitrogen is attached to a long (sphingoid) chain. Here the minimum contiguous chain 
         length (THRESH) is set to 8.
      4) The molecular weight is >500 Da.
      
    Returns:
      bool: True if the molecule is predicted to be a glucosylceramide, False otherwise.
      str: Explanation of the classification decision.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----- 1) Identify exactly one glucose ring -----
    # Use a strict SMARTS including stereochemistry.
    strict_glucose_smarts = "O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    strict_glucose = Chem.MolFromSmarts(strict_glucose_smarts)
    if strict_glucose is None:
        return False, "Error creating strict glucose SMARTS"
    glucose_matches = mol.GetSubstructMatches(strict_glucose)
    
    # If not found, try a looser pattern (ignoring stereochemistry).
    if not glucose_matches:
        loose_glucose_smarts = "OC1OC(O)C(O)C(O)C1O"
        loose_glucose = Chem.MolFromSmarts(loose_glucose_smarts)
        if loose_glucose is None:
            return False, "Error creating loose glucose SMARTS"
        glucose_matches = mol.GetSubstructMatches(loose_glucose)
    
    num_glucose = len(glucose_matches)
    if num_glucose == 0:
        return False, "Glucose head group not found"
    if num_glucose > 1:
        return False, f"Found {num_glucose} sugar rings; expected exactly one glucose head group"
    
    # Record the atom indices of the glucose sugar so that we avoid following into that subgraph.
    sugar_atom_indices = set(glucose_matches[0])
    
    # ----- 2) Detect amide bonds -----
    # SMARTS pattern for a general amide: carbonyl carbon (C=O) adjacent to a nitrogen.
    amide_smarts = "[CX3](=O)[NX3]"
    amide_pattern = Chem.MolFromSmarts(amide_smarts)
    if amide_pattern is None:
        return False, "Error defining amide SMARTS"
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond detected (expected in the ceramide backbone)"
    
    # ----- 3) Helper: DFS to compute longest contiguous chain of carbons (ignoring sugar atoms) -----
    def chain_length(atom, parent_idx, visited, sugar_set):
        """
        Recursively traverse carbon atoms (atomic number 6) ignoring those in sugar_set.
        Only avoid immediately going back to the parent atom.
        """
        best = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == parent_idx:
                continue
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in visited or nbr.GetIdx() in sugar_set:
                continue
            new_visited = visited | {nbr.GetIdx()}
            length = 1 + chain_length(nbr, atom.GetIdx(), new_visited, sugar_set)
            if length > best:
                best = length
        return best

    # ----- 4) Evaluate each amide bond for ceramide-like features -----
    THRESH = 8  # Minimum required contiguous carbon count on each side.
    ceramide_like = False

    # Each amide match is a tuple (carbonyl_idx, nitrogen_idx).
    for match in amide_matches:
        carbonyl_idx, nitrogen_idx = match[0], match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
        
        # For each side, start by “stepping off” from the junction.
        # Acyl (fatty) side: from the carbonyl atom, follow all adjacent carbons except the amide nitrogen.
        max_acyl = 0
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == nitrogen_idx:
                continue  # skip the bond to the amide nitrogen.
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in sugar_atom_indices:
                continue
            visited = {carbonyl_idx, nbr.GetIdx()}
            current = 1 + chain_length(nbr, carbonyl_idx, visited, sugar_atom_indices)
            if current > max_acyl:
                max_acyl = current

        # Sphingoid (base) side: from the nitrogen, follow all adjacent carbons except the carbonyl.
        max_sphingo = 0
        for nbr in nitrogen_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_idx:
                continue
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in sugar_atom_indices:
                continue
            visited = {nitrogen_idx, nbr.GetIdx()}
            current = 1 + chain_length(nbr, nitrogen_idx, visited, sugar_atom_indices)
            if current > max_sphingo:
                max_sphingo = current

        # Debug note: In a true ceramide the fatty acyl chain (from the carbonyl side)
        # and the sphingoid chain (from the nitrogen side) should both be long.
        if max_acyl >= THRESH and max_sphingo >= THRESH:
            ceramide_like = True
            break  # Accept if at least one amide meets the criteria.
    
    if not ceramide_like:
        return False, ("Amide bond not associated with two long aliphatic chains (fatty acyl and sphingoid) "
                       "as expected for a ceramide")
    
    # ----- 5) Final check: Molecular weight should be greater than 500 Da -----
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 500:
        return False, f"Molecular weight too low ({mw:.1f} Da); not typical of a glucosylceramide"
    
    return True, ("Structure has one glucose head group, a ceramide-like amide bond (with long fatty acyl "
                  "and sphingoid chains), and sufficient molecular weight.")

# Example usage (uncomment to test):
# test_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, explanation = is_glucosylceramide(test_smiles)
# print(result, explanation)