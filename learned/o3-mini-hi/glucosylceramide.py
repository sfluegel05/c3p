"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: glucosylceramide
Definition: any cerebroside in which the monosaccharide head group is glucose.
A true glucosylceramide should contain exactly one glucose ring (typically attached via an oxygen),
an amide bond connecting a fatty acyl chain to a sphingoid base, and a high enough molecular weight.
Furthermore, the amide bond should be “ceramide‐like”: the acyl group and the sphingoid part should
each show a contiguous carbon chain of at least a specified length.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    
    Checks:
      1) The molecule contains exactly one glucose ring (first by strict stereochemistry,
         and then by a looser pattern if needed).
      2) The molecule contains at least one amide bond.
      3) At least one amide bond is "ceramide‐like": the carbonyl carbon is attached to one long
         (fatty acyl) chain and the amide nitrogen is attached to another long (sphingoid base) chain.
         Here the minimum contiguous carbon chain length (THRESH) is set to 8.
      4) The molecular weight is >500 Da.
      
    Returns:
      bool: True if the molecule is predicted to be a glucosylceramide, False otherwise.
      str: Explanation of the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----- 1) Find exactly one glucose ring -----
    # Use a strict SMARTS that includes stereochemistry
    strict_glucose_smarts = "O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
    strict_glucose = Chem.MolFromSmarts(strict_glucose_smarts)
    if not strict_glucose:
        return False, "Error creating strict glucose SMARTS"
    glucose_matches = mol.GetSubstructMatches(strict_glucose)
    
    # If strict pattern not found, try a looser one (ignoring stereochemistry)
    if not glucose_matches:
        loose_glucose_smarts = "OC1OC(O)C(O)C(O)C1O"
        loose_glucose = Chem.MolFromSmarts(loose_glucose_smarts)
        if not loose_glucose:
            return False, "Error creating loose glucose SMARTS"
        glucose_matches = mol.GetSubstructMatches(loose_glucose)
    
    num_glucose = len(glucose_matches)
    if num_glucose == 0:
        return False, "Glucose head group not found"
    if num_glucose > 1:
        return False, f"Found {num_glucose} sugar rings; expected exactly one glucose head group"
    
    # Record sugar ring atoms so that we don't follow into that subgraph later.
    sugar_atom_indices = set(glucose_matches[0])
    
    # ----- 2) Check for the presence of an amide bond -----
    # Use a general amide SMARTS pattern
    amide_smarts = "[CX3](=O)[NX3]"
    amide_pattern = Chem.MolFromSmarts(amide_smarts)
    if not amide_pattern:
        return False, "Error defining amide SMARTS"
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond detected (expected in the ceramide backbone)"
    
    # ----- Helper: Compute the longest contiguous carbon chain length -----
    # This DFS traverses only atoms with atomic number == 6 (carbons) and avoids a set of excluded atoms.
    def longest_carbon_chain(atom, visited, exclude_set):
        max_length = 0
        for nbr in atom.GetNeighbors():
            # Only follow carbon atoms that haven't been visited or excluded.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited and nbr.GetIdx() not in exclude_set:
                new_visited = visited | {nbr.GetIdx()}
                branch_length = longest_carbon_chain(nbr, new_visited, exclude_set)
                if branch_length > max_length:
                    max_length = branch_length
        return 1 + max_length  # count current carbon

    # ----- 3) Evaluate each amide bond for ceramide-like features -----
    THRESH = 8  # minimum required contiguous carbons on each side
    ceramide_like = False

    # For each matched amide bond, interpret match[0] as the carbonyl carbon and match[1] as the amide nitrogen.
    for match in amide_matches:
        carbonyl_idx, nitrogen_idx = match[0], match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
        # For chain search avoid the sugar atoms and the amide atoms themselves.
        avoid_set = sugar_atom_indices | {carbonyl_idx, nitrogen_idx}
        
        # --- Fatty acyl side: from the carbonyl carbon, select neighbors other than nitrogen and oxygen (i.e. not the C=O oxygen)
        acyl_lengths = []
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == nitrogen_idx:
                continue  # skip the amide nitrogen
            # Skip heteroatoms (e.g. carbonyl oxygen)
            if nbr.GetAtomicNum() != 6:
                continue
            length = longest_carbon_chain(nbr, {nbr.GetIdx()}, avoid_set)
            acyl_lengths.append(length)
        max_acyl = max(acyl_lengths) if acyl_lengths else 0
        
        # --- Sphingoid side: from the amide nitrogen, select neighbors other than the carbonyl carbon.
        sphingo_lengths = []
        for nbr in nitrogen_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_idx:
                continue  # do not go back to the carbonyl
            if nbr.GetAtomicNum() != 6:
                continue
            length = longest_carbon_chain(nbr, {nbr.GetIdx()}, avoid_set)
            sphingo_lengths.append(length)
        max_sphingo = max(sphingo_lengths) if sphingo_lengths else 0
        
        # Check whether both sides have long enough chains.
        if max_acyl >= THRESH and max_sphingo >= THRESH:
            ceramide_like = True
            break

    if not ceramide_like:
        return False, ("Amide bond not associated with two long aliphatic chains (fatty acyl and sphingoid) "
                       "as expected for a ceramide")

    # ----- 4) Final check: molecular weight should be greater than 500 Da -----
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 500:
        return False, f"Molecular weight too low ({mw:.1f} Da); not typical of a glucosylceramide"
    
    return True, "Structure has one glucose head group, a ceramide-like amide (with long acyl and sphingoid chains), and sufficient molecular weight."

# Example usage:
# test_smiles = "CCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, explanation = is_glucosylceramide(test_smiles)
# print(result, explanation)