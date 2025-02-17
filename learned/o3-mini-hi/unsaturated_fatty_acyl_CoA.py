"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: unsaturated fatty acyl-CoA

Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of an unsaturated fatty acid.
Improvement over previous version:
 - Enforces the presence of an adenine ring (common to CoA) to better recognize the full CoA moiety.
 - Requires that the fatty acyl fragment (detached from the thioester group) contains at least 6 carbon atoms.
 - Inspects the fatty acyl fragment for at least one non‐aromatic C=C bond.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is an unsaturated fatty acyl-CoA.
    
    Criteria:
       1. The molecule must contain a CoA-like moiety. In addition to a partial pattern like 
          "SCCNC(=O)CCNC(=O)" seen in many acyl-CoAs, we now require the adenine ring fragment.
       2. The molecule must contain a thioester group linking an acyl fragment to the CoA portion.
       3. The fatty acyl chain (the side attached to the thioester carbonyl via its R group)
          must be “long” (≥6 carbon atoms) and contain at least one non‐aromatic C=C double bond.
          
    Args:
       smiles (str): SMILES string for the molecule.
       
    Returns:
       (bool, str): (True, reason) if the molecule is classified as an unsaturated fatty acyl-CoA,
                    (False, reason) otherwise.
    """
    # Parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Check for CoA moiety ---
    # Old check was based on a fragment "SCCNC(=O)CCNC(=O)".
    # We now also require an adenine fragment, common to CoA (pattern taken from many CoA structures).
    coa_core_smarts = "SCCNC(=O)CCNC(=O)"
    adenine_smarts = "n1cnc2ncnc12"  # a pattern for the adenine ring
    coa_core = Chem.MolFromSmarts(coa_core_smarts)
    adenine_core = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(coa_core) or not mol.HasSubstructMatch(adenine_core):
        return False, "Coenzyme A moiety not detected"
    
    # --- Search for thioester group ---
    # Look for a carbonyl (C=O) directly attached to a sulfur.
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found linking acyl chain to CoA"
    
    # --- For each thioester group, analyze the acyl chain side ---
    # In a thioester group R-C(=O)S-, the fatty acyl chain is the R group attached to the carbonyl.
    for match in thioester_matches:
        # match is a tuple of atom indices corresponding roughly to (carbonyl carbon, sulfur)
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Find the neighbor on the carbonyl which is not the oxygen (double-bonded) and not the sulfur.
        fatty_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue
            # Skip the oxygen (part of the C=O)
            if nbr.GetAtomicNum() == 8:
                continue
            # Expect a carbon (for fatty acyl chain)
            if nbr.GetAtomicNum() == 6:
                fatty_neighbor = nbr.GetIdx()
                break
        
        if fatty_neighbor is None:
            continue  # This thioester does not lead to a plausible fatty acyl chain
        
        # To focus on the fatty acyl fragment, break the bond from carbonyl to the acyl chain.
        rwmol = Chem.RWMol(mol)
        try:
            rwmol.RemoveBond(carbonyl_idx, fatty_neighbor)
        except Exception:
            continue  # if removal fails, move to next thioester candidate
        
        # Identify fragments (each fragment is given as a tuple of atom indices)
        frags = rdmolops.GetMolFrags(rwmol, asMols=False)
        fatty_frag = None
        for frag in frags:
            if fatty_neighbor in frag:
                fatty_frag = frag
                break
        if fatty_frag is None:
            continue
        
        # Count carbon atoms in the fatty acyl fragment.
        nC = sum(1 for idx in fatty_frag if rwmol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if nC < 6:
            # Too short to be a fatty acyl chain
            continue
        
        # Check for unsaturation (at least one non‐aromatic C=C double bond) in the fatty fragment.
        unsaturation_found = False
        for bond in rwmol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            # Consider bond only if both atoms are part of the fatty fragment.
            if a1 in fatty_frag and a2 in fatty_frag:
                if (bond.GetBondType() == Chem.BondType.DOUBLE and
                    not bond.GetIsAromatic()):
                    # Ensure both atoms are carbons (expected for simple acyl chains)
                    atom1 = rwmol.GetAtomWithIdx(a1)
                    atom2 = rwmol.GetAtomWithIdx(a2)
                    if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                        unsaturation_found = True
                        break
        
        if unsaturation_found:
            return True, "Contains CoA moiety with thioester-linked unsaturated fatty acyl chain"
    
    # If none of the thioester groups provided a fatty fragment with unsaturation:
    return False, "Thioester group found but fatty acyl chain is either too short or appears fully saturated"

# For testing/debugging purposes, you might run the function if this file is executed directly.
if __name__ == '__main__':
    # Example: (5E)-tetradecenoyl-CoA from the provided data
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)CCC/C=C/CCCCCCCC)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_unsaturated_fatty_acyl_CoA(test_smiles)
    print(result, reason)