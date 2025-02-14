"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A BCFA is characterized by the presence of a carboxylic acid group (or its ionized form),
    and one or more branching points typically involving smaller alkyl groups attached 
    to the main hydrocarbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: A tuple containing a boolean value indicating if it's a BCFA, and the reason
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for carboxylic acid group (including ionized form)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O,H]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Find branching points (carbon attached to at least 3 other carbons)
    branch_patterns = [Chem.MolFromSmarts("[C;D3]"), Chem.MolFromSmarts("[C;D4]")]
    branch_found = False
    for pattern in branch_patterns:
        for match in mol.GetSubstructMatches(pattern):
            atom = mol.GetAtomWithIdx(match[0])
            if atom.GetDegree() > 2:  # Implies branching
                if not any(atom.IsInRingSize(size) for size in (3,)):  # Cyclopropyl ring exception
                    branch_found = True
                    break
    
    if not branch_found:
        return False, "No suitable branching pattern found"

    return True, "Contains a branched hydrocarbon chain with carboxylic acid group"

# Test examples (provide one or two example SMILES strings here to check classification)
smiles_examples = [
    "CC(C)CCCCCCC(=O)O",  # Expected: True, branched chain
    "CCCCCCCCCCCCCCCCCC(=O)O",  # Expected: False, straight chain
]

# Output results for testing examples
for sm in smiles_examples:
    result, reason = is_branched_chain_fatty_acid(sm)
    print(f"SMILES: {sm} -> Result: {result}, Reason: {reason}")