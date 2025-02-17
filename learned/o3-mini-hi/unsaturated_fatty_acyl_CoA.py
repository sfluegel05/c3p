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
  - Breaking the bond to isolate the fatty acyl fragment (excluding the CoA part) using fragmentation.
  - Using a SMARTS query to detect a non‐aromatic C=C double bond and counting carbon atoms in the fragment.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is an unsaturated fatty acyl-CoA.
    
    The molecule is classified as an unsaturated fatty acyl-CoA if it meets these criteria:
       1. Contains an adenine ring (indicating the CoA moiety).
       2. Contains a thioester linkage ([CX3](=O)[SX2]).
       3. The acyl fragment (obtained by disconnecting the thioester bond) is sufficiently long 
          (at least 6 carbon atoms) and has at least one non‐aromatic C=C double bond.
          
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       (bool, str): (True, reason) if the molecule is an unsaturated fatty acyl-CoA,
                    (False, reason) otherwise.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Check for CoA moiety via adenine substructure ---
    # A very simple adenine SMARTS pattern.
    adenine_smarts = "c1nc2c(n1)nc(nc2)"
    adenine_query = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(adenine_query):
        return False, "Coenzyme A moiety not detected (adenine ring missing)"
    
    # --- Look for thioester linkage: carbonyl directly bound to a sulfur ---
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester_query = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_query)
    if not thioester_matches:
        return False, "No thioester group found linking acyl chain to CoA"
    
    # Prepare a SMARTS query to detect non‐aromatic C=C bond in a fatty chain
    unsat_query = Chem.MolFromSmarts("C=C")
    
    # Process each thioester match to try to extract a fatty acyl fragment
    for match in thioester_matches:
        # match is a tuple: (carbonyl_idx, sulfur_idx)
        carbonyl_idx, sulfur_idx = match[0], match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the fatty branch: a neighbor of the carbonyl that is not the sulfur
        fatty_neighbor_idx = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue  # skip sulfur
            # Also skip oxygen atoms (e.g. carbonyl oxygen)
            if nbr.GetAtomicNum() == 8:
                continue
            # We expect a carbon that will lead to a fatty chain
            if nbr.GetAtomicNum() == 6:
                fatty_neighbor_idx = nbr.GetIdx()
                break
        if fatty_neighbor_idx is None:
            continue  # try next thioester candidate
        
        # Attempt to break the bond between the carbonyl and the fatty branch.
        try:
            # Note: attach dummy atoms (*) at the cleavage sites.
            fragmented = Chem.FragmentOnBonds(mol, [(carbonyl_idx, fatty_neighbor_idx)], addDummies=True)
        except Exception as e:
            continue
        
        # Get the disconnected fragments as separate molecules.
        frags = rdmolops.GetMolFrags(fragmented, asMols=True, sanitizeFrags=True)
        if not frags or len(frags) < 2:
            continue
        
        # Select the fragment that does NOT contain the adenine substructure.
        fatty_candidates = []
        for frag in frags:
            if not frag.HasSubstructMatch(adenine_query):
                fatty_candidates.append(frag)
        if not fatty_candidates:
            continue
        
        # From the candidates, choose the fragment with the largest number of carbon atoms.
        best_frag = None
        max_carbons = 0
        for frag in fatty_candidates:
            # Count carbon atoms.
            c_count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
            if c_count > max_carbons:
                max_carbons = c_count
                best_frag = frag
        
        # Check that the chain is long enough.
        if max_carbons < 6:
            # Not a sufficiently long fatty acyl chain.
            continue
        
        # Use a SMARTS query to check if the best fragment contains a non‐aromatic C=C double bond.
        if best_frag.HasSubstructMatch(unsat_query):
            return True, ("Contains adenine moiety linked via a thioester to a fatty acyl chain "
                          "with at least 6 carbons and a non‐aromatic C=C double bond")
    
    # If none of the thioester candidates produced a valid fatty acyl fragment:
    return False, "Thioester group found but fatty acyl chain is either too short or fully saturated"

# If executed as a script, run an example test.
if __name__ == '__main__':
    # Test one of the provided examples: (5E)-tetradecenoyl-CoA
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)CCC/C=C/CCCCCCCC)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_unsaturated_fatty_acyl_CoA(test_smiles)
    print(result, reason)