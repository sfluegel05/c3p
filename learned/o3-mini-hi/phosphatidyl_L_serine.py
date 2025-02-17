"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
#!/usr/bin/env python
"""
Classifies: Phosphatidyl-L-serine (an aminophospholipid in which a phosphatidyl group is esterified to the hydroxy group of serine).
We now use a more specific substructure SMARTS that embeds the key connectivity:
    P(=O)(O)(O)OC[C](N)C(=O)O
This should only match molecules in which a phosphate group is directly attached via an oxygen to the CH2 of serine.
"""

from rdkit import Chem

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    
    A phosphatidyl-L-serine must have a phosphoserine headgroup. In phosphatidyl-L-serine the
    serine headgroup (HOCH2CH(N)COOH) becomes esterified on its hydroxyl (CH2OH becomes O–)
    by a phosphatidyl group. Thus we look for a substructure in which a phosphate group (P with one =O and two –OH/O–)
    is attached via an oxygen (O) to a CH2 group that is then connected to a NH and C(=O)O fragment.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phosphatidyl-L-serine, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the input SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the phosphoserine headgroup.
    # This pattern encodes a phosphate group with one double-bonded oxygen and two other oxygen substituents,
    # where one oxygen is esterified to a CH2 group that bears a serine unit (CH2CH(N)C(=O)O).
    # We do not enforce chirality here for flexibility.
    phosphoserine_smarts = "P(=O)(O)(O)OC[C](N)C(=O)O"
    ps_pattern = Chem.MolFromSmarts(phosphoserine_smarts)
    if ps_pattern is None:
        return False, "Error creating SMARTS pattern for phosphoserine"
    
    # Check if the molecule contains the phosphoserine headgroup
    if mol.HasSubstructMatch(ps_pattern):
        return True, "Found phosphoserine headgroup: phosphate esterified to serine"
    else:
        return False, "No phosphoserine headgroup (P(=O)(O)(O)OC[C](N)C(=O)O) found"

# (Optional) Testing examples
if __name__ == "__main__":
    # Example: 1-oleoyl-2-arachidonoyl-sn-glycero-3-phospho-L-serine
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC"
    result, reason = is_phosphatidyl_L_serine(test_smiles)
    print(result, reason)