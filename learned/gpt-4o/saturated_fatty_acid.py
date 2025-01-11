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
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid functional group found"
    
    # Check for any unsaturated bonds (C=C or C#C)
    unsaturated_double_pattern = Chem.MolFromSmarts("[C]=[C]")
    unsaturated_triple_pattern = Chem.MolFromSmarts("[C]#[C]")
    if mol.HasSubstructMatch(unsaturated_double_pattern) or mol.HasSubstructMatch(unsaturated_triple_pattern):
        return False, "Contains carbon-carbon multiple bonds (unsaturated)"

    # Acknowledge that ring structures may be present, but are less common
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 1:
        return False, "Contains complex ring structures; simplest saturated fatty acids do not"

    # Check for sufficient length of carbon chain to be considered a fatty acid
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, "Carbon chain too short for a fatty acid"
    
    # It passed all tests to be considered a saturated fatty acid
    return True, "Contains carboxylic acid group, is linear/branched, and contains no unsaturations"

# Test the function with examples
examples = [
    "CCCCC(O)=O",  # valeric acid
    "CCC(C)CCCCCCCCCCCCCCCCC(O)=O",  # 18-methylicosanoic acid
    "CCCCCCCCCCCCCCCCCCC(O)=O",  # nonadecanoic acid
    "CC(C)CCCCCC(O)=O",  # 7-methyloctanoic acid
    "CCCC(CCC)C(O)=O",  # valproic acid
]

for example in examples:
    result, reason = is_saturated_fatty_acid(example)
    print(f"SMILES: {example}, Result: {result}, Reason: {reason}")