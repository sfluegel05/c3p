"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Cyclic Fatty Acid
Definition: Any fatty acid containing anywhere in its structure a ring of atoms.
A fatty acid is expected to have a carboxylic acid group and a long aliphatic (or unsaturated) chain.
We add a heuristic that the molecule must have a contiguous chain of 5 or more non‐ring, non‐aromatic carbon atoms.
"""

from rdkit import Chem

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    
    A cyclic fatty acid is defined as a fatty acid having a carboxyl group,
    a long (at least 5-carbon) acyclic chain, and at least one ring.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a cyclic fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a broader SMARTS pattern for the carboxylic acid group.
    # This pattern will match both protonated [OX2H] and deprotonated [OX1] forms.
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
        
    # Check if the molecule contains any rings.
    ring_info = mol.GetRingInfo()
    if not ring_info or len(ring_info.AtomRings()) == 0:
        return False, "No ring found in the structure"
        
    # Use a SMARTS pattern to detect a contiguous chain of at least 5 non‐ring, non‐aromatic carbons.
    # This is a heuristic for a "fatty" chain.
    chain_pattern = Chem.MolFromSmarts("[#6;!R;!a]{5,}")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No sufficiently long acyclic aliphatic chain found"

    # Passed all checks: carboxylic acid, ring, and long chain.
    return True, "Contains carboxylic acid group, a long aliphatic chain, and at least one ring in the structure"


# For debugging and testing purposes:
if __name__ == "__main__":
    # (R)-lipoic acid: known cyclic fatty acid example.
    test_smiles = "OC(=O)CCCC[C@@H]1CCSS1"
    result, reason = is_cyclic_fatty_acid(test_smiles)
    print(result, reason)