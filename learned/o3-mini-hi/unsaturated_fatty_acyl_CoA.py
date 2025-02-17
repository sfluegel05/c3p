"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: unsaturated fatty acyl-CoA

Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of an unsaturated fatty acid.

Improvement over previous version:
 - Instead of requiring an exact match to a short CoA fragment, we now demand a proper adenine
   ring (a core component of CoA) to be present. This increases the chance to detect the full CoA moiety.
 - We require that a thioester group linking the acyl chain to CoA is present.
 - We ensure the detached acyl chain has at least 6 carbons and bears at least one non‐aromatic C=C bond.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is an unsaturated fatty acyl-CoA.

    Criteria:
       1. The molecule must contain a CoA-like moiety. Here we require at least the presence of an adenine ring.
       2. The molecule must contain a thioester group linking an acyl fragment to the CoA portion.
       3. The acyl chain attached to the thioester carbonyl (R group) must have at least 6 carbons
          and contain at least one non‐aromatic C=C double bond.
    
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
    # Instead of using a rigid fragment, we require that the molecule contains an adenine ring,
    # which is an essential part of the CoA structure.
    # Here we use a SMARTS pattern for adenine.
    adenine_smarts = "c1nc2c(n1)nc(nc2)"
    adenine_core = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(adenine_core):
        return False, "Coenzyme A moiety not detected (adenine ring missing)"
    
    # --- Search for thioester group ---
    # Look for a carbonyl (C=O) directly attached to a sulfur, i.e. the thioester linkage.
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found linking acyl chain to CoA"
    
    # --- Analyze each thioester group candidate ---
    for match in thioester_matches:
        # match is a tuple of atom indices corresponding roughly to (carbonyl carbon, sulfur)
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the neighbor on the carbonyl that is not the doubly bonded oxygen and not the sulfur.
        fatty_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == sulfur_idx:
                continue
            # Skip any oxygen (assumed to form the carbonyl group)
            if nbr.GetAtomicNum() == 8:
                continue
            if nbr.GetAtomicNum() == 6:
                fatty_neighbor = nbr.GetIdx()
                break
        if fatty_neighbor is None:
            continue  # no plausible fatty acyl side detected for this thioester
        
        # Copy the molecule into an editable RWMol and remove the bond between the carbonyl and the fatty fragment.
        rwmol = Chem.RWMol(mol)
        try:
            rwmol.RemoveBond(carbonyl_idx, fatty_neighbor)
        except Exception:
            continue  # if bond removal fails, try next thioester group
        
        # Identify fragments: one fragment should be the detached acyl chain.
        frags = rdmolops.GetMolFrags(rwmol, asMols=False)
        fatty_frag = None
        for frag in frags:
            if fatty_neighbor in frag:
                fatty_frag = frag
                break
        if fatty_frag is None:
            continue
        
        # Check that the fatty acyl fragment is "long" enough (at least 6 carbons)
        n_carbons = sum(1 for idx in fatty_frag if rwmol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if n_carbons < 6:
            continue
        
        # Check for at least one non‐aromatic carbon–carbon double bond in the fatty acid fragment.
        unsaturation_found = False
        for bond in rwmol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            # Consider bond only if both atoms are part of the fatty fragment.
            if a1 in fatty_frag and a2 in fatty_frag:
                if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic():
                    atom1 = rwmol.GetAtomWithIdx(a1)
                    atom2 = rwmol.GetAtomWithIdx(a2)
                    # Check that both atoms are carbons
                    if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                        unsaturation_found = True
                        break
        
        if unsaturation_found:
            return True, "Contains adenine (CoA) moiety with thioester-linked unsaturated fatty acyl chain"
            
    # If no valid fatty acyl fragment was found among the thioester candidates:
    return False, "Thioester group found but fatty acyl chain is either too short or fully saturated"

# For testing/debugging purposes, the function below may be executed directly.
if __name__ == '__main__':
    # Example: (5E)-tetradecenoyl-CoA (first provided example)
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)CCC/C=C/CCCCCCCC)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_unsaturated_fatty_acyl_CoA(test_smiles)
    print(result, reason)