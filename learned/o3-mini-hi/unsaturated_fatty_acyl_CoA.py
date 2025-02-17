"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: unsaturated fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any unsaturated fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is an unsaturated fatty acyl-CoA.
    
    Criteria:
      1. The molecule must contain a CoA moiety. We check for a characteristic fragment 
         (e.g. "SCCNC(=O)CCNC(=O)") seen in many CoA derivatives.
      2. The molecule must contain a thioester group ([CX3](=O)[SX2]) linking the fatty acyl chain to CoA.
      3. The fatty acyl chain (the part attached to the thioester carbonyl carbon, on the “R” side)
         must contain at least one non‐aromatic carbon-carbon double bond, in order to be “unsaturated”.
         
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): True with a success message if the molecule is classified as an unsaturated 
                   fatty acyl-CoA; False with the reason otherwise.
    """
    # Parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, check for a CoA fragment. We use a characteristic SMARTS common in CoA:
    # (This pattern is not perfect but appears in many acyl-CoA structures.)
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not detected"
    
    # Next, search for a thioester group: a carbonyl carbon attached to a sulfur.
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found linking acyl chain to CoA"
    
    # For each thioester group found, check whether its acyl side is unsaturated.
    # In a thioester group R-C(=O)S-, the fatty acyl chain is the R attached to the carbonyl.
    for match in thioester_matches:
        # match is a tuple of atom indices corresponding to [C(=O)S]
        carbonyl_idx = match[0]
        sulfur_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the neighbor attached to the carbonyl that is not the oxygen (C=O)
        # nor the sulfur (which leads to the CoA portion)
        fatty_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            # Skip if the neighbor is the sulfur (part of thioester)
            if nbr.GetIdx() == sulfur_idx:
                continue
            # Skip the oxygen (the carbonyl oxygen)
            if nbr.GetAtomicNum() == 8:
                continue
            # We prefer a carbon neighbor (expected for a fatty acyl chain)
            if nbr.GetAtomicNum() == 6:
                fatty_neighbor = nbr.GetIdx()
                break
        
        if fatty_neighbor is None:
            continue  # try the next thioester, if any
        
        # Now, to study only the fatty acyl chain, we remove (temporarily) the bond from the
        # carbonyl to the fatty acyl chain. We make a copy of the molecule for fragmentation.
        rwmol = Chem.RWMol(mol)
        try:
            rwmol.RemoveBond(carbonyl_idx, fatty_neighbor)
        except Exception:
            continue  # if bond removal fails, try next match
        
        # Identify fragments. GetMolFrags returns tuples of atom indices.
        frags = rdmolops.GetMolFrags(rwmol, asMols=False)
        # Find the fragment that contains the fatty_neighbor index.
        fatty_frag = None
        for frag in frags:
            if fatty_neighbor in frag:
                fatty_frag = frag
                break
        
        if fatty_frag is None:
            continue
        
        # Now, inspect all bonds in the fatty acyl fragment.
        # We want to find at least one non-aromatic C=C double bond.
        unsaturation_found = False
        for bond in rwmol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            # Consider only bonds where both atoms are in the fatty fragment.
            if a1 in fatty_frag and a2 in fatty_frag:
                # Check if the bond is a double bond, non-aromatic, and connects two carbons.
                if (bond.GetBondType() == Chem.BondType.DOUBLE and 
                    not bond.GetIsAromatic()):
                    atom1 = rwmol.GetAtomWithIdx(a1)
                    atom2 = rwmol.GetAtomWithIdx(a2)
                    if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                        unsaturation_found = True
                        break
        
        if unsaturation_found:
            return True, "Contains CoA moiety with thioester-linked unsaturated fatty acyl chain"
    
    # If none of the thioester groups yield an unsaturated fatty acyl chain,
    return False, "Thioester group found, but fatty acyl chain appears fully saturated"

# For testing/debugging purposes:
if __name__ == '__main__':
    # Example SMILES for (5E)-tetradecenoyl-CoA (one of the provided examples)
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)CCC/C=C/CCCCCCCC)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    result, reason = is_unsaturated_fatty_acyl_CoA(test_smiles)
    print(result, reason)