"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: short-chain fatty acyl-CoA

A short-chain fatty acyl-CoA is defined as a fatty acyl-CoA whose acyl portion 
(derived from a short-chain fatty acid, typically 2–6 carbons in length) is attached 
via a thioester bond (C(=O)–S) to a Coenzyme A (CoA) moiety.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA.
    
    The classification checks for:
      - A unique thioester bond (SMARTS: [CX3](=O)[SX2]).
      - The presence of a CoA moiety (simplified by matching a characteristic substructure).
      - That after breaking the thioester bond, the acyl fragment (identified via a marked 
        carbonyl atom) contains between 2 and 6 carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a short-chain fatty acyl-CoA, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES into an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS for the thioester bond: carbonyl (C=O) linked to sulfur (S)
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester_pat = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_pat)
    if not thioester_matches:
        return False, "No thioester bond (C(=O)-S) found"
    if len(thioester_matches) > 1:
        return False, f"Multiple thioester bonds found ({len(thioester_matches)}); expected exactly one"

    # In the match tuple, assume first atom is the carbonyl carbon and second is the sulfur.
    carbonyl_idx, sulfur_idx = thioester_matches[0][0], thioester_matches[0][1]

    # Check for the presence of a CoA moiety using a simplified SMARTS pattern.
    # The chosen pattern is a characteristic fragment that should be part of CoA.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pat = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pat):
        return False, "Coenzyme A moiety not detected"

    # Mark the carbonyl atom with a property so we can track which fragment is the acyl chain.
    atom = mol.GetAtomWithIdx(carbonyl_idx)
    atom.SetProp("is_acyl_marker", "true")

    # Create an editable copy of the molecule.
    rw_mol = Chem.RWMol(mol)
    # Verify that the bond between the carbonyl and sulfur exists.
    bond = rw_mol.GetBondBetweenAtoms(carbonyl_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond not found for removal"
    # Remove the thioester bond to disconnect the acyl chain from the CoA moiety.
    rw_mol.RemoveBond(carbonyl_idx, sulfur_idx)

    # Get the disconnected fragments as separate molecules.
    frags = Chem.GetMolFrags(rw_mol, asMols=True, sanitizeFrags=True)
    if len(frags) != 2:
        return False, f"Expected 2 fragments (acyl and CoA), but got {len(frags)} fragments"

    # Identify the acyl fragment as the fragment that contains the marked carbonyl atom.
    acyl_frag = None
    coa_frag = None
    for frag in frags:
        for atom in frag.GetAtoms():
            if atom.HasProp("is_acyl_marker"):
                acyl_frag = frag
                break
        # If not acyl, check for CoA marker by substructure match.
        if acyl_frag is None and frag.HasSubstructMatch(coa_pat):
            coa_frag = frag

    if acyl_frag is None:
        # In case marker was not located (could be due to fragmentation issues)
        return False, "Acyl fragment not identified via marked carbonyl"
    
    # Count carbon atoms in the acyl fragment.
    acyl_carbon_count = sum(1 for atom in acyl_frag.GetAtoms() if atom.GetAtomicNum() == 6)
    if acyl_carbon_count < 2 or acyl_carbon_count > 6:
        return False, f"Acyl chain has {acyl_carbon_count} carbon(s); expected a short-chain fatty acyl group (2–6 carbons)"
    
    return True, f"Contains a thioester bond with a short-chain fatty acyl group ({acyl_carbon_count} carbons) attached to a CoA moiety"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Testing one of the examples: 2-methylbutanoyl-CoA
    test_smiles = "CCC(C)C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    result, reason = is_short_chain_fatty_acyl_CoA(test_smiles)
    print("Result:", result)
    print("Reason:", reason)