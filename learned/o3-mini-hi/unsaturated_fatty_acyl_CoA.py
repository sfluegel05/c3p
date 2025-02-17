"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: unsaturated fatty acyl-CoA

Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of an unsaturated fatty acid.

Improvements over previous version:
 - Instead of matching a rigid CoA fragment, we now demand an adenine substructure.
 - We search for the thioester linkage ([CX3](=O)[SX2]) and then break the bond between the carbonyl carbon and its fatty (acyl) side.
 - We classify the fatty acyl fragment by requiring that it be “long” (at least 6 carbons) and contain at least one non‐aromatic C=C double bond.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is an unsaturated fatty acyl-CoA.
    
    The molecule is classified as an unsaturated fatty acyl-CoA if it meets these criteria:
       1. Contains an adenine ring (CoA moiety indicator).
       2. Contains a thioester linkage [CX3](=O)[SX2].
       3. The acyl fragment attached to the carbonyl (i.e. the fatty chain) contains at least 6 carbons 
          and at least one non-aromatic C=C double bond.
    
    Args:
       smiles (str): SMILES string for the molecule.
       
    Returns:
       (bool, str): (True, reason) if the molecule is classified as an unsaturated fatty acyl-CoA,
                    (False, reason) otherwise.
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Check for a CoA moiety via the adenine substructure ---
    adenine_smarts = "c1nc2c(n1)nc(nc2)"  # a minimal adenine pattern
    adenine_query = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(adenine_query):
        return False, "Coenzyme A moiety not detected (adenine ring missing)"
    
    # --- Look for thioester linkage: carbonyl directly bound to a sulfur ---
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester_query = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_query)
    if not thioester_matches:
        return False, "No thioester group found linking acyl chain to CoA"
    
    # For each thioester candidate, try to extract the fatty acyl fragment.
    for match in thioester_matches:
        # match gives indices: (carbonyl carbon, sulfur)
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the neighbor of the carbonyl that is not the sulfur (this should be the fatty acyl branch).
        fatty_neighbor_idx = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue
            # Skip oxygens (part of the carbonyl group).
            if nbr.GetAtomicNum() == 8:
                continue
            if nbr.GetAtomicNum() == 6:
                fatty_neighbor_idx = nbr.GetIdx()
                break
        if fatty_neighbor_idx is None:
            continue  # no plausible fatty acyl branch found
        
        # Use FragmentOnBonds to break the bond between the carbonyl and the fatty branch.
        try:
            fragmented = Chem.FragmentOnBonds(mol, [ (carbonyl_idx, fatty_neighbor_idx) ], addDummies=True)
        except Exception:
            continue
        frags = rdmolops.GetMolFrags(fragmented, asMols=True, sanitizeFrags=True)
        if not frags or len(frags) < 2:
            continue
        
        # We assume one fragment is the CoA moiety (contains adenine) and the other is the fatty chain.
        fatty_frag = None
        for frag in frags:
            if not frag.HasSubstructMatch(adenine_query):
                fatty_frag = frag
                break
        
        if fatty_frag is None:
            continue  # could not separate a fatty chain from the CoA moiety.
        
        # Count the number of carbon atoms in the fatty fragment.
        n_carbons = sum(1 for atom in fatty_frag.GetAtoms() if atom.GetAtomicNum() == 6)
        if n_carbons < 6:
            continue  # fatty chain too short
        
        # Check for at least one non‐aromatic C=C double bond within the fatty fragment.
        unsaturation_found = False
        for bond in fatty_frag.GetBonds():
            # Only consider bonds between two carbon atoms.
            if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic():
                a1 = bond.GetBeginAtom()
                a2 = bond.GetEndAtom()
                if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                    unsaturation_found = True
                    break
        
        if unsaturation_found:
            return True, "Contains adenine (CoA) moiety with thioester-linked unsaturated fatty acyl chain"
    
    # If no valid fatty acyl fragment is identified:
    return False, "Thioester group found but fatty acyl chain is either too short or fully saturated"


# For testing/debugging purposes the function below may be executed directly.
if __name__ == '__main__':
    # Test examples provided (example: (5E)-tetradecenoyl-CoA)
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)CCC/C=C/CCCCCCCC)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_unsaturated_fatty_acyl_CoA(test_smiles)
    print(result, reason)