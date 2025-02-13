"""
Classifies: CHEBI:75874 para-terphenyl
"""
from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule belongs to the class of para-terphenyl compounds based on its SMILES string.
    A para-terphenyl is characterized by a structure consisting of three phenyl rings connected in a linear fashion (para orientation).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a para-terphenyl or its derivative, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for para-terphenyl
    # The pattern represents three connected phenyl rings with potential substitutions
    para_terphenyl_pattern = Chem.MolFromSmarts('c1cc(ccc1)-c2ccc(cc2)-c3ccc(cc3)')

    if not mol.HasSubstructMatch(para_terphenyl_pattern):
        return False, "Does not match para-terphenyl core structure"

    return True, "Contains para-terphenyl core structure"

# Example test
test_smiles = "c1cc(ccc1)-c2ccc(cc2)-c3ccc(cc3)"  # Simplified para-terphenyl structure
result, reason = is_para_terphenyl(test_smiles)
print(f"Result: {result}, Reason: {reason}")