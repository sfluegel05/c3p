"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: Very long‐chain fatty acyl‐CoA
Definition: A fatty acyl‐CoA in which the fatty acyl group (the acyl chain attached via a thioester bond)
has a chain length greater than C22 (i.e. 23 or more carbon atoms in the acyl fragment).
Our strategy:
  1. Check for the thioester group using C(=O)[S] SMARTS.
  2. Check for a CoA moiety by looking for an adenine fragment.
  3. Mark the carbonyl carbon, break (fragment) the molecule at the bond to sulfur and
     then count the number of carbon atoms in the fragment that contains that marked atom.
  4. If the count is at least 23, then we classify the acyl-CoA as very long-chain.
"""

from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    Criteria:
      1. The molecule must have a thioester group (a carbonyl bonded to a sulfur atom).
      2. The molecule must contain a CoA moiety (detected via adenine substructure patterns).
      3. After “cutting” the molecule at the thioester bond (between the carbonyl carbon and sulfur),
         the fatty acyl chain (i.e. the fragment containing the carbonyl carbon) must contain at least 23 carbons.
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool, str: (True, explanation) if the molecule meets the criteria,
                   (False, explanation) otherwise.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # STEP 1: Look for a thioester group.
    # Pattern: a carbonyl carbon (C=O) directly bonded to a sulfur atom.
    thioester_smarts = "C(=O)[S]"  
    thioester_pat = Chem.MolFromSmarts(thioester_smarts)
    if thioester_pat is None:
        return False, "Error creating thioester SMARTS pattern"
    thioester_matches = mol.GetSubstructMatches(thioester_pat)
    if not thioester_matches:
        return False, "No thioester group found; not an acyl-CoA"
    
    # For simplicity, use the first found match.
    # According to our SMARTS, match[0] is the carbonyl carbon and match[1] is the sulfur.
    carbonyl_idx, sulfur_idx = thioester_matches[0][0], thioester_matches[0][1]
    
    # STEP 2: Check for CoA moiety.
    # We look for adenine fragments which are common in CoA.
    adenine_smarts1 = "n1cnc2c(N)ncnc12"
    adenine_smarts2 = "n1cnc2ncnc12"
    adenine_pat1 = Chem.MolFromSmarts(adenine_smarts1)
    adenine_pat2 = Chem.MolFromSmarts(adenine_smarts2)
    if adenine_pat1 is None or adenine_pat2 is None:
        return False, "Error creating adenine SMARTS patterns"
    if not (mol.HasSubstructMatch(adenine_pat1) or mol.HasSubstructMatch(adenine_pat2)):
        return False, "No CoA moiety detected (adenine fragment missing)"
    
    # STEP 3: Mark the carbonyl atom and break the thioester bond.
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    # Mark it with a property so we can trace it after fragmentation.
    carbonyl_atom.SetProp("is_acyl", "1")
    
    # Find the bond between the carbonyl carbon and the sulfur.
    bond = mol.GetBondBetweenAtoms(carbonyl_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond not found between carbonyl and sulfur"
    bond_idx = bond.GetIdx()

    # We'll fragment the molecule by breaking the thioester bond.
    # The resulting molecule will have dummy atoms (atomic number 0) at the break points.
    try:
        frag_mol = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
    except Exception as e:
        return False, f"Error fragmenting molecule: {str(e)}"
    
    # Get the fragments (as separate molecules).
    frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments obtained after bond breaking"
    
    # Locate the fragment that contains the acyl chain.
    # It must have the marked carbonyl atom (which we set the property "is_acyl").
    acyl_frag = None
    for frag in frags:
        for atom in frag.GetAtoms():
            if atom.HasProp("is_acyl") and atom.GetProp("is_acyl") == "1":
                acyl_frag = frag
                break
        if acyl_frag is not None:
            break
    
    if acyl_frag is None:
        return False, "Could not isolate fatty acyl fragment from the thioester bond"

    # STEP 4: Count the number of carbon atoms in the acyl fragment.
    # Dummy atoms (atomic number 0) will be skipped.
    carbon_count = sum(1 for atom in acyl_frag.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # For a fatty acyl chain having a chain length > C22, we require at least 23 carbon atoms.
    if carbon_count < 23:
        return False, f"Fatty acyl chain has {carbon_count} carbons; must be at least 23"
    
    return True, f"Fatty acyl chain has {carbon_count} carbons, which is >=23"

# (Optional) Example usage:
# test_smiles = ("CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)"
#                "[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)"
#                "[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
# result, reason = is_very_long_chain_fatty_acyl_CoA(test_smiles)
# print(result, reason)