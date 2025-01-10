"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid contains a carboxylic acid group and no carbon-carbon multiple bonds,
    typically has a linear or branched alkyl chain structure, and typically contains a moderate to long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,H0-]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid functional group found"
    
    # Check for any unsaturated bonds (C=C or C#C)
    unsaturated_pattern = Chem.MolFromSmarts("[C]=[C] | [C]#[C]")
    if mol.HasSubstructMatch(unsaturated_pattern):
        return False, "Contains carbon-carbon multiple bonds (unsaturated)"

    # Check for presence of cycle which is uncommon in saturated fatty acids
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains ring structure; typical saturated fatty acids do not"

    # Check for sufficient length of carbon chain to be considered a fatty acid
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, "Carbon chain too short for a fatty acid"
    
    # It passed all tests to be considered a saturated fatty acid
    return True, "Contains carboxylic acid group, is linear/branched, and contains no unsaturations"


# Test the function with provided examples
examples = [
    "CCCCC(O)=O",  # valeric acid
    "CCC(C)CCCCCCCCCCCCCCCCC(O)=O",  # 18-methylicosanoic acid
    "CCCCCCCCCCCCCCCCCCC(O)=O",  # nonadecanoic acid
    "CC(C)CCCCCC(O)=O",  # 7-methyloctanoic acid
    # Add additional examples if needed
]

for example in examples:
    result, reason = is_saturated_fatty_acid(example)
    print(f"SMILES: {example}, Result: {result}, Reason: {reason}")