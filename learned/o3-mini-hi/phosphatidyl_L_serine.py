"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
#!/usr/bin/env python
"""
Classifies: Phosphatidyl-L-serine
A phosphatidyl-L-serine is an aminophospholipid in which a phosphatidyl group is esterified to the hydroxy group of serine.
Thus it contains a phosphoserine headgroup (i.e. a phosphate attached via an oxygen to a CH2 that is further linked to NH and C(=O)O) along with two fatty acyl (ester) chains.
This classifier uses not only a SMARTS match for the key phosphoserine headgroup but also enforces that
• the molecule has a molecular weight high enough to support two acyl chains,
• and that at least two ester bond fragments (typical of acyl chains) are found.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    
    The algorithm checks for:
      1. A phosphoserine headgroup. We use the SMARTS "P(=O)(O)(O)OC[C](N)C(=O)O" pattern which encodes:
         a phosphate group (one double-bonded O and two –OH/O– substituents) attached to a CH2,
         which is then connected to a serine fragment (CH(N)C(=O)O).
      2. That the molecule has a high molecular weight (>600 Da), typical for the full phospholipid that contains 
         two fatty acyl (lipid) chains.
      3. At least two long acyl chains detected via the presence of ester bonds connected to a long alkyl chain.
         We search using a SMARTS pattern that looks for an ester function "OC(=O)CCCC", which requires 
         an ester bond followed by at least four carbon atoms (this is a proxy for a fatty acyl chain).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a phosphatidyl-L-serine, False otherwise.
        str: Explanation of the classification decision.
    """
    
    # Parse the input SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for the phosphoserine headgroup.
    # This pattern encodes: phosphate (P(=O)(O)(O)) attached via O to CH2, then CH(N) and C(=O)O.
    phosphoserine_smarts = "P(=O)(O)(O)OC[C](N)C(=O)O"
    ps_pattern = Chem.MolFromSmarts(phosphoserine_smarts)
    if ps_pattern is None:
        return False, "Error creating SMARTS pattern for phosphoserine"
    
    # Check for the phosphoserine headgroup
    if not mol.HasSubstructMatch(ps_pattern):
        return False, "No phosphoserine headgroup (P(=O)(O)(O)OC[C](N)C(=O)O) found"
    
    # Check the molecular weight.
    # Full phosphatidyl-L-serines (with two fatty acyl chains) are typically heavy (>600 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for full phosphatidyl-L-serine"
    
    # Define a SMARTS pattern for an acyl chain ester bond.
    # We look for an ester bond "OC(=O)" that is attached to an alkyl fragment of at least 4 carbons (CCCC).
    # This is used as a proxy for a fatty acyl chain.
    acyl_smarts = "OC(=O)CCCC"
    acyl_pattern = Chem.MolFromSmarts(acyl_smarts)
    if acyl_pattern is None:
        return False, "Error creating SMARTS pattern for acyl chains"
    
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    
    # We require at least two matches to ensure that two fatty acyl chains are present.
    if len(acyl_matches) < 2:
        return False, f"Only found {len(acyl_matches)} acyl chain ester(s); expected at least 2 for phosphatidyl-L-serine"
    
    # If all checks pass, classify the molecule as a phosphatidyl-L-serine.
    return True, "Found phosphoserine headgroup with sufficient molecular weight and acyl chains"

# (Optional) Testing examples when run as a script.
if __name__ == "__main__":
    # Example: 1-oleoyl-2-arachidonoyl-sn-glycero-3-phospho-L-serine
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC"
    result, reason = is_phosphatidyl_L_serine(test_smiles)
    print(result, reason)