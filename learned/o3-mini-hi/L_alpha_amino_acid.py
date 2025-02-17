"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
#!/usr/bin/env python
"""
Classifies: L-alpha-amino acid
Definition: Any alpha-amino acid having L-configuration at the alpha-carbon.
An L-alpha-amino acid should contain an amino acid backbone N–C(α)(sidechain)–C(=O)O
with the chiral center (C(α)) having L configuration.
"""

from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    
    An L-alpha-amino acid must have an amino acid backbone, i.e., an alpha-carbon 
    bonded to a nitrogen and a carboxylate (C(=O)O) group, and the alpha-carbon must be chiral 
    with the proper configuration corresponding (by CIP rules) to the L form (typically 'S').
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise.
        str : Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemical information is computed
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define SMARTS patterns for the alpha-amino acid backbone.
    # This pattern looks for an amino group bonded to a chiral carbon that is in turn bonded to a carboxylic acid.
    # We prepare two patterns to cover both possible chirality flag styles in the SMILES:
    pattern1 = Chem.MolFromSmarts("N[C@H]([*])(C(=O)O)")
    pattern2 = Chem.MolFromSmarts("N[C@@H]([*])(C(=O)O)")
    
    # Try to find matches using one of the patterns.
    matches = mol.GetSubstructMatches(pattern1) + mol.GetSubstructMatches(pattern2)
    if not matches:
        return False, "Alpha-amino acid backbone (N-C(alpha)-C(=O)O) not found"
    
    # Check each match to verify the alpha-carbon has an assigned CIP code that corresponds to L configuration
    # (by convention for most amino acids under CIP, L is 'S').
    for match in matches:
        # match is a tuple of atom indices from the SMARTS;
        # we expect: index0: amino nitrogen, index1: alpha-carbon, index2: carboxyl carbon.
        alpha_idx = match[1]
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        # Check that the alpha carbon is chiral and has a computed CIP code.
        if not alpha_atom.HasProp("_CIPCode"):
            continue  # no CIP code assigned so cannot determine configuration
        cip_code = alpha_atom.GetProp("_CIPCode")
        if cip_code == "S":
            return True, "Found an alpha-amino acid backbone with L (S CIP) configuration at the alpha-carbon"
    
    return False, "Alpha-amino acid backbone found but alpha-carbon does not have L configuration (S CIP code)"
    
# Example usage (uncomment for testing):
# test_smiles = "N[C@@H](CC(=O)c1cc(O)ccc1N)C(O)=O"  # 5-hydroxy-L-kynurenine example
# result, reason = is_L_alpha_amino_acid(test_smiles)
# print(result, reason)