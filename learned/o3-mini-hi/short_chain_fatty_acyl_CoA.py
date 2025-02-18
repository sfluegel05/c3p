"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: short-chain fatty acyl-CoA

A short-chain fatty acyl-CoA is defined as a fatty acyl-CoA whose acyl portion (derived from a short-chain fatty acid, typically 2–6 carbons in length) is attached via a thioester bond (C(=O)–S) to a coenzyme A (CoA) moiety.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA.
    
    The classification checks for:
      - A thioester bond that indicates an acyl binding (pattern: [CX3](=O)[SX2])
      - The presence of a CoA moiety (recognized via a substructure match, here simplified as "SCCNC(=O)CCNC(=O)")
      - An acyl chain (the fragment which results when breaking the bond between the acyl carbon and sulfur)
        having 2 to 6 carbon atoms (as is common for short-chain fatty acids).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a short-chain fatty acyl-CoA, False otherwise.
        str: Reason for the classification decision.
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a thioester bond: carbonyl (C=O) linked to sulfur (S)
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester_pat = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pat)
    if not thioester_matches:
        return False, "No thioester bond (C(=O)–S) found"
    if len(thioester_matches) > 1:
        return False, f"Multiple thioester bonds found ({len(thioester_matches)}); expected exactly one for a fatty acyl-CoA"
    
    # Get the unique thioester match.
    # The SMARTS "[CX3](=O)[SX2]" should return a tuple (C_atom, S_atom),
    # where the C atom is the one in the carbonyl group.
    carbonyl_idx, sulfur_idx = thioester_matches[0][0], thioester_matches[0][1]

    # Check for the presence of a CoA moiety.
    # We use a simplified substructure pattern typical for CoA fragments found in our examples.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pat = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pat):
        return False, "Coenzyme A moiety not detected"

    # To isolate the acyl group, break the bond between the carbonyl carbon and sulfur.
    # The acyl group is the fragment that contains the carbonyl carbon after removing the thioester bond.
    rw_mol = Chem.RWMol(mol)
    rw_mol.RemoveBond(carbonyl_idx, sulfur_idx)

    # Obtain the fragments (each fragment is returned as a tuple of original atom indices)
    frags = Chem.GetMolFrags(rw_mol, asMols=False)
    acyl_frag_indices = None
    for frag in frags:
        if carbonyl_idx in frag:
            acyl_frag_indices = frag
            break
    if acyl_frag_indices is None:
        return False, "Failed to isolate acyl fragment"

    # Count the number of carbon atoms in the acyl fragment.
    # This count (including the carbonyl carbon) should be 2 to 6 for a short-chain fatty acid.
    acyl_carbon_count = 0
    for idx in acyl_frag_indices:
        atom = rw_mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:  # Carbon
            acyl_carbon_count += 1

    if acyl_carbon_count < 2 or acyl_carbon_count > 6:
        return False, f"Acyl chain has {acyl_carbon_count} carbon(s); expected a short-chain fatty acyl group (2–6 carbons)"

    return True, f"Contains a thioester bond with a short-chain fatty acyl group ({acyl_carbon_count} carbons) attached to a CoA moiety"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test with one of the provided examples: isovaleryl-CoA
    test_smiles = "CC(C)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_short_chain_fatty_acyl_CoA(test_smiles)
    print("Result:", result)
    print("Reason:", reason)