"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: CHEBI short-chain fatty acyl-CoA

Definition:
  A short-chain fatty acyl-CoA is defined as a fatty acyl-CoA where the acyl part 
  (coming from a short-chain fatty acid, typically 2–6 carbons including the carbonyl carbon)
  is attached via a thioester bond (C(=O)-S) to a coenzyme A (CoA) moiety.
  
Algorithm:
  1. Parse the SMILES.
  2. Search for a unique thioester bond using the SMARTS pattern "[CX3](=O)[SX2]".
  3. Use FragmentOnBonds (with addDummies=True) to “cut” the molecule at the thioester bond.
  4. Identify the acyl fragment by counting its carbon atoms (ignoring dummy atoms).
     The acyl chain (including its carbonyl carbon) should contain 2–6 C atoms.
"""

from rdkit import Chem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    
    The classification checks for:
      - A unique thioester bond (SMARTS: [CX3](=O)[SX2]).
      - Fragmentation at that bond, and that the fragment corresponding to the acyl group
        contains between 2 and 6 carbon atoms (including the carbonyl carbon).
      - Implicitly, the remainder should represent a coenzyme A moiety.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a short-chain fatty acyl-CoA, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the thioester bond (C(=O)-S)
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester_pat = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pat)
    if not thioester_matches:
        return False, "No thioester bond (C(=O)-S) found"
    if len(thioester_matches) > 1:
        return False, f"Multiple thioester bonds found ({len(thioester_matches)}); expected exactly one"
    
    # Get the atoms from the unique match: assume first atom is the carbonyl carbon and second is sulfur.
    carbonyl_idx, sulfur_idx = thioester_matches[0][0], thioester_matches[0][1]
    bond = mol.GetBondBetweenAtoms(carbonyl_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond not found in the molecule"
    
    # Now fragment the molecule at the thioester bond. We use addDummies=True so that dummy atoms [*]
    # mark the cut sites and the fragments remain intact.
    frag_mol = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=True)
    frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
    
    if len(frags) != 2:
        return False, f"Expected 2 fragments after bond cleavage, but got {len(frags)}"
    
    # For each fragment, count the number of carbon atoms; ignore dummy atoms (atomic number 0).
    def count_carbons(mol_fragment):
        return sum(1 for atom in mol_fragment.GetAtoms() if atom.GetAtomicNum() == 6)
    
    frag_carbon_counts = [count_carbons(frag) for frag in frags]
    
    # The acyl chain (including its carbonyl carbon) should have between 2 and 6 carbons.
    acyl_frag = None
    acyl_count = None
    for frag, count in zip(frags, frag_carbon_counts):
        # If the fragment has a low carbon count, consider it the acyl chain.
        if 2 <= count <= 6:
            acyl_frag = frag
            acyl_count = count
            break

    if acyl_frag is None:
        return False, f"Acyl chain not identified; fragment carbon counts found: {frag_carbon_counts}"

    # For a bit more clarity, provide the fragment counts.
    return True, f"Contains a unique thioester bond. The acyl chain fragment has {acyl_count} carbon(s) (expected 2–6)."

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test with one of the provided examples: 2-methylbutanoyl-CoA
    test_smiles = "CCC(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_short_chain_fatty_acyl_CoA(test_smiles)
    print("Result:", result)
    print("Reason:", reason)