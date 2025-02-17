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
    
    # attempt to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----- Check for a single glucose ring head group -----
    # First try a strict SMARTS which encodes stereochemistry
    strict_glucose_smarts = "O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O"
    strict_glucose_pattern = Chem.MolFromSmarts(strict_glucose_smarts)
    if strict_glucose_pattern is None:
        return False, "Error in defining strict glucose SMARTS"
    glucose_matches = mol.GetSubstructMatches(strict_glucose_pattern)
    
    # If that fails, try a looser pattern ignoring stereochemistry
    if not glucose_matches:
        loose_glucose_smarts = "OC1OC(O)C(O)C(O)C1O"
        loose_glucose_pattern = Chem.MolFromSmarts(loose_glucose_smarts)
        if loose_glucose_pattern is None:
            return False, "Error in defining loose glucose SMARTS"
        glucose_matches = mol.GetSubstructMatches(loose_glucose_pattern)
    num_glucose = len(glucose_matches)
    if num_glucose == 0:
        return False, "Glucose head group not found"
    if num_glucose > 1:
        return False, f"Found {num_glucose} sugar rings; expected exactly one glucose head group"
    
    # Record the atoms that form the sugar
    sugar_atom_indices = set( glucose_matches[0] )
    
    # ----- Check for an amide bond -----
    # We use a general amide SMARTS
    amide_smarts = "[CX3](=O)[NX3]"
    amide_pattern = Chem.MolFromSmarts(amide_smarts)
    if amide_pattern is None:
        return False, "Error in defining amide SMARTS pattern"
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond detected (expected in the ceramide backbone)"
    
    # Helper: recursively compute the maximum length of a contiguous chain of carbon atoms.
    # We count only carbon atoms (atomic number 6) and avoid atoms that are in 'exclude_set'.
    def dfs_chain(atom, visited, exclude_set):
        best = 0
        for nbr in atom.GetNeighbors():
            # only count carbon neighbors not already visited and not in the exclude set
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited and nbr.GetIdx() not in exclude_set:
                new_visited = visited | {nbr.GetIdx()}
                length = dfs_chain(nbr, new_visited, exclude_set)
                if length > best:
                    best = length
        return 1 + best  # count the current carbon as well

    # We now check each amide match to see if it looks like part of a ceramide backbone.
    # For each match the pattern gives a tuple (carbonyl carbon index, amide nitrogen index).
    # On the carbonyl side we expect a fatty acyl chain, and on the N side a sphingoid chain.
    # We require each chain to have a maximum contiguous chain length of at least THRESH (here 8).
    THRESH = 8
    ceramide_like = False
    for match in amide_matches:
        carbonyl_idx, nitrogen_idx = match[0], match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
        
        # From the carbonyl atom, get the neighbor (other than the amide N and excluding oxygen) 
        # that is likely the start of an acyl (fatty acid) chain.
        acyl_candidate = None
        for nbr in carbonyl_atom.GetNeighbors():
            # Skip the amide nitrogen and the oxygen (the carbonyl oxygen) 
            if nbr.GetIdx() == nitrogen_idx:
                continue
            if nbr.GetAtomicNum() == 8:
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_candidate = nbr
                break
        if acyl_candidate is None:
            continue  # no candidate acyl chain
        
        # For the nitrogen, get the neighbor (other than the carbonyl carbon) that is carbon.
        sphingo_candidate = None
        for nbr in nitrogen_atom.GetNeighbors():
            if nbr.GetIdx() == carbonyl_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                sphingo_candidate = nbr
                break
        if sphingo_candidate is None:
            continue  # no sphingoid chain candidate
        
        # Compute longest contiguous chain lengths from the candidate atoms;
        # we force the DFS to ignore sugar atoms (and also the amide atoms themselves to avoid backtracking into the acyl head)
        acyl_chain_length = dfs_chain(acyl_candidate, {acyl_candidate.GetIdx()}, sugar_atom_indices | {carbonyl_idx, nitrogen_idx})
        sphingo_chain_length = dfs_chain(sphingo_candidate, {sphingo_candidate.GetIdx()}, sugar_atom_indices | {carbonyl_idx, nitrogen_idx})
        
        if acyl_chain_length >= THRESH and sphingo_chain_length >= THRESH:
            ceramide_like = True
            break  # found at least one ceramide-like amide bond

    if not ceramide_like:
        return False, "Amide bond not associated with two long aliphatic chains (fatty acyl and sphingoid) as expected for a ceramide"
    
    # ----- Check molecular weight (typical glucosylceramides are large) -----
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 500:
        return False, f"Molecular weight too low ({mw:.1f} Da); not typical of a glucosylceramide"
    
    return True, "Structure has a single glucose head group, a ceramide-like amide (with long acyl and sphingoid chains), and sufficient molecular weight."

# Example: (uncomment below lines to test with one SMILES example)
# example_smiles = "CCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_glucosylceramide(example_smiles)
# print(result, reason)