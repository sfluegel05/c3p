"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid contains a carboxylic acid group and no carbon-carbon multiple bonds.

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

    # Check for the presence of carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid functional group found"
        
    # Check for any carbon-carbon double bonds (C=C)
    carbon_double_bond_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(carbon_double_bond_pattern):
        return False, "Contains carbon-carbon double bonds (unsaturated)"

    # Check for any carbon-carbon triple bonds (C#C)
    carbon_triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if mol.HasSubstructMatch(carbon_triple_bond_pattern):
        return False, "Contains carbon-carbon triple bonds (unsaturated)"

    # Ensure there is a significant carbon chain (at least 4 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Carbon chain too short for fatty acid"

    return True, "Contains carboxylic acid group and no carbon-carbon multiple bonds"


# You can test the function with the given examples
examples = [
    "CCCCC(O)=O",  # valeric acid
    "CCC(C)CCCCCCCCCCCCCCCCC(O)=O",  # 18-methylicosanoic acid
    # Add more examples if needed
]

for example in examples:
    result, reason = is_saturated_fatty_acid(example)
    print(f"SMILES: {example}, Result: {result}, Reason: {reason}")