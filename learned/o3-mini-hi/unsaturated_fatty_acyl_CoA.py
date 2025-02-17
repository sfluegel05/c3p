"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: unsaturated fatty acyl-CoA

Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any unsaturated fatty acid.

This version improves upon the previous code by:
  - Checking for the adenine substructure as an indicator of CoA.
  - Locating the thioester linkage using SMARTS.
  - Breaking the bond between the carbonyl and fatty acyl branch and selecting the fragment 
    (excluding the CoA moiety) that has the largest number of carbon atoms.
  - Ensuring that the fatty chain is sufficiently long (>=6 carbons) and contains at least one 
    non‐aromatic C=C double bond.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is an unsaturated fatty acyl-CoA.
    
    The molecule is classified as an unsaturated fatty acyl-CoA if it meets these criteria:
       1. Contains an adenine ring (indicating the CoA moiety).
       2. Contains a thioester linkage ([CX3](=O)[SX2]).
       3. The acyl fragment attached to the carbonyl contains at least 6 carbons 
          and at least one non-aromatic C=C double bond.
          
    Returns:
       (bool, str): (True, reason) if the molecule is an unsaturated fatty acyl-CoA,
                    (False, reason) otherwise.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Check for a CoA moiety via the adenine substructure ---
    # This minimal adenine SMARTS helps indicate the CoA part.
    adenine_smarts = "c1nc2c(n1)nc(nc2)"
    adenine_query = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(adenine_query):
        return False, "Coenzyme A moiety not detected (adenine ring missing)"
    
    # --- Look for thioester linkage: carbonyl directly bound to a sulfur ---
    # SMARTS for thioester: a carbon (sp2) double bonded to O and single bonded to S.
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester_query = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_query)
    if not thioester_matches:
        return False, "No thioester group found linking acyl chain to CoA"
    
    # We will try each thioester found to see if we can extract a qualifying fatty acyl fragment.
    for match in thioester_matches:
        # match returns a tuple: (carbonyl carbon, sulfur)
        carbonyl_idx, sulfur_idx = match[0], match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify fatty acyl branch: a neighbor of the carbonyl that is not the sulfur (and not oxygen of C=O)
        fatty_neighbor_idx = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue
            if nbr.GetAtomicNum() == 8:  # skip oxygen (e.g. carbonyl O)
                continue
            if nbr.GetAtomicNum() == 6:  # carbon candidate
                fatty_neighbor_idx = nbr.GetIdx()
                break
        if fatty_neighbor_idx is None:
            continue  # no plausible fatty branch found here
        
        # Attempt to fragment (break bond) between the carbonyl and its fatty neighbor.
        try:
            fragmented = Chem.FragmentOnBonds(mol, [(carbonyl_idx, fatty_neighbor_idx)], addDummies=True)
        except Exception:
            continue
        # Get the disconnected fragments (as separate molecule objects).
        frags = rdmolops.GetMolFrags(fragmented, asMols=True, sanitizeFrags=True)
        if not frags or len(frags) < 2:
            continue
        
        # From the fragments, select those that do not contain adenine (i.e. likely the fatty chain).
        fatty_candidates = []
        for frag in frags:
            if not frag.HasSubstructMatch(adenine_query):
                fatty_candidates.append(frag)
        if not fatty_candidates:
            continue

        # Among candidate fragments, choose the one with the highest carbon count.
        best_frag = None
        max_carbons = 0
        for frag in fatty_candidates:
            n_carbons = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
            if n_carbons > max_carbons:
                max_carbons = n_carbons
                best_frag = frag

        # Ensure the fatty chain is long enough.
        if max_carbons < 6:
            # Fatty candidate too short
            continue

        # Check for at least one non‐aromatic C=C double bond in the fatty fragment.
        unsaturation_found = False
        for bond in best_frag.GetBonds():
            # Check only bonds between two carbons.
            if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic():
                a1 = bond.GetBeginAtom()
                a2 = bond.GetEndAtom()
                if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                    unsaturation_found = True
                    break
        
        if unsaturation_found:
            return True, "Contains adenine moiety linked via a thioester to a sufficiently long unsaturated fatty acyl chain"
    
    # If we evaluated all thioester candidates without success:
    return False, "Thioester group found but fatty acyl chain is either too short or fully saturated"


# If executed as a script, run an example test.
if __name__ == '__main__':
    # Test one of the provided examples: (5E)-tetradecenoyl-CoA
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)CCC/C=C/CCCCCCCC)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_unsaturated_fatty_acyl_CoA(test_smiles)
    print(result, reason)