"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
#!/usr/bin/env python
"""
Classifies: Phosphatidyl-L-serine (an aminophospholipid where a phosphatidyl group is esterified to the hydroxy group of serine)
"""

from rdkit import Chem

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    
    A phosphatidyl-L-serine must have a serine headgroup with the pattern –OCC(N)C(=O)O
    and the oxygen (that would be esterified in this motif) must be bonded to a phosphorus.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phosphatidyl-L-serine, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS for the serine head group (ignoring explicit chirality)
    # This pattern matches –OCC(N)C(=O)O which is typical for the serine headgroup in phosphatidyl-L-serine.
    serine_smarts = "[OX2]CC(N)C(=O)[O]"
    serine_pattern = Chem.MolFromSmarts(serine_smarts)
    
    # Check if the serine motif exists in the molecule.
    matches = mol.GetSubstructMatches(serine_pattern)
    if not matches:
        return False, "No serine headgroup motif (-OCC(N)C(=O)O) found"
    
    # For each match, check if the oxygen (the first atom in the match) is directly bonded to a phosphorus.
    for match in matches:
        # match[0] corresponds to the oxygen (the -O that is the start of the serine group)
        oxygen_atom = mol.GetAtomWithIdx(match[0])
        # Loop over the neighbors of the oxygen to see if any is P (atomic number 15)
        for neighbor in oxygen_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 15:
                return True, "Found serine head group with oxygen esterified to a phosphate (P)"
    
    # If no oxygen of the serine groups is directly attached to phosphorus, then it is not a phosphatidyl species.
    return False, "Serine headgroup present but not esterified to a phosphate"

# (Optional) For testing purposes one could call the function here:
if __name__ == "__main__":
    # Example SMILES for 1-oleoyl-2-arachidonoyl-sn-glycero-3-phospho-L-serine
    test_smiles = "CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC"
    result, reason = is_phosphatidyl_L_serine(test_smiles)
    print(result, reason)