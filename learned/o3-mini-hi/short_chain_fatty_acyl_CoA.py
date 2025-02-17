"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: CHEBI short-chain fatty acyl-CoA

Definition:
  A short-chain fatty acyl-CoA is defined as a fatty acyl-CoA where the acyl part 
  (coming from a short-chain fatty acid, typically 2–6 carbons) is attached via a 
  thioester bond (i.e. C(=O)-S) to a coenzyme A (CoA) moiety.

The algorithm implemented here:
  1. Parses the SMILES.
  2. Searches for a unique thioester bond defined by SMARTS "[CX3](=O)[SX2]".
  3. Checks for the presence of a characteristic CoA fragment via a simplified SMARTS pattern.
  4. Uses RDKit’s FragmentOnBonds (with addDummies=False) to "cut" the thioester bond.
  5. Obtains and assigns the resulting fragments as the acyl part (short, 2–6 carbons)
     or the CoA moiety.
"""

from rdkit import Chem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    
    The classification checks for:
      - A unique thioester bond (SMARTS: [CX3](=O)[SX2]).
      - The presence of a CoA moiety (via a characteristic substructure search).
      - That cleavage of the thioester bond produces one fragment (the acyl chain)
        with between 2 and 6 carbon atoms.
    
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
    
    # Obtain the bond index for the thioester bond
    # We assume that the first two atoms in the match correspond to 
    # the carbonyl carbon and the sulfur.
    carbonyl_idx, sulfur_idx = thioester_matches[0][0], thioester_matches[0][1]
    bond = mol.GetBondBetweenAtoms(carbonyl_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond not found in the molecule"
    
    # Check for the presence of a CoA fragment via a simplified pattern.
    # (This pattern is not exhaustive but expected to occur in the CoA portion.)
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pat = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pat):
        return False, "Coenzyme A moiety not detected"
    
    # Use FragmentOnBonds to cleave the molecule at the thioester bond.
    # Using addDummies=False yields chemically sensible fragments.
    frag_mol = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=False)
    frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
    
    if len(frags) != 2:
        return False, f"Expected 2 fragments after bond cleavage, but got {len(frags)}"
    
    # Identify which fragment is the CoA moiety by checking for the CoA substructure,
    # and which fragment is the acyl chain.
    coa_frag = None
    acyl_frag = None
    for frag in frags:
        if frag.HasSubstructMatch(coa_pat):
            coa_frag = frag
        else:
            acyl_frag = frag

    # If we did not get a clear assignment, then try to use the number of carbons:
    if acyl_frag is None:
        for frag in frags:
            carbon_count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
            if 2 <= carbon_count <= 6:
                acyl_frag = frag
            else:
                coa_frag = frag

    # Safety check for acyl fragment identification.
    if acyl_frag is None:
        counts = [sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6) for frag in frags]
        return False, f"Acyl chain not identified; fragment carbon counts found: {counts}"
    
    # Count the number of carbon atoms (atomic number 6) in the acyl fragment.
    acyl_carbon_count = sum(1 for atom in acyl_frag.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (2 <= acyl_carbon_count <= 6):
        return False, f"Acyl chain has {acyl_carbon_count} carbon(s); expected a short-chain fatty acyl group (2–6 carbons)"
    
    return True, f"Contains a thioester bond with a short-chain fatty acyl group ({acyl_carbon_count} carbons) attached to a CoA moiety"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Using one of the provided examples: 2-methylbutanoyl-CoA
    test_smiles = "CCC(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_short_chain_fatty_acyl_CoA(test_smiles)
    print("Result:", result)
    print("Reason:", reason)