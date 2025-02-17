"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: glucosylceramide
Definition: any cerebroside in which the monosaccharide head group is glucose.
A true glucosylceramide should contain exactly one glucose ring (typically attached via an oxygen),
an amide bond connecting a fatty acyl chain to a sphingoid base, and a sufficiently high molecular weight.
Furthermore, the amide bond should be “ceramide‐like”: the acyl group and the sphingoid part should appear as long aliphatic chains.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on a SMILES string.
    
    Checks:
      1) The molecule contains exactly one glucose ring. First using a strict sugar SMARTS
         (which includes stereochemistry), and if that fails then a looser SMARTS.
      2) The molecule contains an amide bond [CX3](=O)[NX3].
      3) One of those amide bonds appears to be part of a ceramide: namely, the carbonyl carbon must be
         connected to a sufficiently long aliphatic “fatty acid” chain and the amide nitrogen must be attached
         to a sufficiently long sphingoid base chain. (Here we require a chain of at least 8 contiguous carbons.)
      4) The molecular weight is >500 Da.
      
    Returns:
      bool: True if the molecule is predicted to be a glucosylceramide, False otherwise.
      str: Explanation of the classification decision.
    """
    
    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----- Check for the presence of exactly one glucose ring -----
    # First, try a strict SMARTS (includes stereochemistry)
    strict_glucose_smarts = "O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O"
    strict_glucose_pattern = Chem.MolFromSmarts(strict_glucose_smarts)
    if strict_glucose_pattern is None:
        return False, "Error defining strict glucose SMARTS"
    glucose_matches = mol.GetSubstructMatches(strict_glucose_pattern)
    
    # If no match is found, try a looser SMARTS (ignoring stereochemistry)
    if not glucose_matches:
        loose_glucose_smarts = "OC1OC(O)C(O)C(O)C1O"
        loose_glucose_pattern = Chem.MolFromSmarts(loose_glucose_smarts)
        if loose_glucose_pattern is None:
            return False, "Error defining loose glucose SMARTS"
        glucose_matches = mol.GetSubstructMatches(loose_glucose_pattern)
    
    num_glucose = len(glucose_matches)
    if num_glucose == 0:
        return False, "Glucose head group not found"
    if num_glucose > 1:
        return False, f"Found {num_glucose} sugar rings; expected exactly one glucose head group"
    
    # Record the atom indices that form the sugar to avoid following into the sugar region later.
    sugar_atom_indices = set(glucose_matches[0])
    
    # ----- Check for an amide bond (general pattern) -----
    amide_smarts = "[CX3](=O)[NX3]"
    amide_pattern = Chem.MolFromSmarts(amide_smarts)
    if amide_pattern is None:
        return False, "Error defining amide SMARTS"
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond detected (expected in the ceramide backbone)"
    
    # Helper function: recursively compute maximum contiguous chain length along carbon atoms,
    # avoiding atoms in exclude_set and already visited atoms.
    def dfs_chain(atom, visited, exclude_set):
        max_length = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited and nbr.GetIdx() not in exclude_set:
                new_visited = visited | {nbr.GetIdx()}
                branch_length = dfs_chain(nbr, new_visited, exclude_set)
                if branch_length > max_length:
                    max_length = branch_length
        return 1 + max_length  # count the current carbon

    # ----- Evaluate each amide bond for a ceramide-like feature -----
    THRESH = 8  # minimum contiguous carbon chain length required on both sides
    ceramide_like = False
    
    for match in amide_matches:
        # In the SMARTS, the first atom is the carbonyl carbon and the second is the amide nitrogen.
        carbonyl_idx, nitrogen_idx = match[0], match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
        
        # For the acyl (fatty acid) side: look at carbonyl carbon neighbors other than:
        #   (a) the amide nitrogen, and (b) any oxygen (i.e. the carbonyl oxygen)
        acyl_chain_lengths = []
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == nitrogen_idx:
                continue
            if nbr.GetAtomicNum() == 8:  # likely the carbonyl oxygen
                continue
            if nbr.GetAtomicNum() == 6:
                length = dfs_chain(nbr, {nbr.GetIdx()}, sugar_atom_indices | {carbonyl_idx, nitrogen_idx})
                acyl_chain_lengths.append(length)
        max_acyl_length = max(acyl_chain_lengths) if acyl_chain_lengths else 0
        
        # For the sphingoid side: from the nitrogen, get carbon neighbors excluding the carbonyl carbon.
        sphingo_chain_lengths = []
        for nbr in nitrogen_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                length = dfs_chain(nbr, {nbr.GetIdx()}, sugar_atom_indices | {carbonyl_idx, nitrogen_idx})
                sphingo_chain_lengths.append(length)
        max_sphingo_length = max(sphingo_chain_lengths) if sphingo_chain_lengths else 0
        
        if max_acyl_length >= THRESH and max_sphingo_length >= THRESH:
            ceramide_like = True
            break  # found a ceramide-like amide bond

    if not ceramide_like:
        return False, ("Amide bond not associated with two long aliphatic chains (fatty acyl and sphingoid) "
                       "as expected for a ceramide")
    
    # ----- Final check: molecular weight should be >500 Da -----
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 500:
        return False, f"Molecular weight too low ({mw:.1f} Da); not typical of a glucosylceramide"
    
    return True, "Structure has a single glucose head group, a ceramide-like amide (with long acyl and sphingoid chains), and sufficient molecular weight."

# Example usage (uncomment the following lines for testing):
# example_smiles = "CCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_glucosylceramide(example_smiles)
# print(result, reason)