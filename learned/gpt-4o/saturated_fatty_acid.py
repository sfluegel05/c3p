"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid contains a carboxylic acid group and no carbon-carbon multiple bonds,
    typically exhibits a linear or branched alkyl chain structure.

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

    # Check for the presence of carboxylic acid group (-C(=O)[O;H1,H0-])
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,H0-]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid functional group found"

    # Check for any carbon-carbon double bonds (C=C) or triple bonds (C#C)
    unsaturated_pattern = Chem.MolFromSmarts("C=C | C#C")
    if mol.HasSubstructMatch(unsaturated_pattern):
        return False, "Contains carbon-carbon multiple bonds (unsaturated)"

    # Check for linear or branched alkyl chain pattern for the backbone
    alkyl_chain_pattern = Chem.MolFromSmarts("C(-C)(-C)(-C)") # Represents a simple branching chain or linear chain
    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "Does not resemble a typical fatty acid alkyl chain"

    # While initially considering carbon count for fatty acids, exceptions exist with shorter carbon count
    # Decide based on literary definitions and the provided known examples which might include smaller acids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short and not a known exception for a saturated fatty acid"
    
    return True, "Contains carboxylic acid group and resembles a typical fatty acid alkyl chain"


# You can test the function with the given examples
examples = [
    "CCCCC(O)=O",  # valeric acid
    "CCC(C)CCCCCCCCCCCCCCCCC(O)=O",  # 18-methylicosanoic acid
    # Add more examples if needed
]

for example in examples:
    result, reason = is_saturated_fatty_acid(example)
    print(f"SMILES: {example}, Result: {result}, Reason: {reason}")