"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Cyclic Fatty Acid
Definition: Any fatty acid containing anywhere in its structure a ring of atoms.
A fatty acid should contain a carboxylic acid group (C(=O)O) and, 
for cyclic fatty acids, at least one ring.
"""

from rdkit import Chem

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    
    A cyclic fatty acid is defined as a fatty acid (containing a carboxyl group)
    that also contains one or more rings anywhere in its structure.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a cyclic fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a SMARTS pattern for the carboxylic acid group.
    # This pattern looks for a carbon double-bonded to oxygen and single-bonded to an -OH.
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
        
    # Check if the molecule contains any rings.
    ring_info = mol.GetRingInfo()
    if not ring_info or len(ring_info.AtomRings()) == 0:
        return False, "No ring found in the structure"
        
    # If the molecule passes both tests, then it is a cyclic fatty acid.
    return True, "Contains carboxylic acid group and at least one ring in the structure"
    
# Example usage (for debugging purposes):
if __name__ == "__main__":
    test_smiles = "OC(=O)CCCC[C@@H]1CCSS1"  # (R)-lipoic acid, a known cyclic fatty acid
    result, reason = is_cyclic_fatty_acid(test_smiles)
    print(result, reason)