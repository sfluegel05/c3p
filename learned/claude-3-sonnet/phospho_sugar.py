"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: CHEBI:27838 phospho sugar
A phospho sugar is any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phosphate group patterns
    phosphate_patterns = [Chem.MolFromSmarts("[P+](O)(O)(O)=O"),
                          Chem.MolFromSmarts("OP(O)(O)=O"),
                          Chem.MolFromSmarts("P(=O)(O)(O)"),
                          Chem.MolFromSmarts("P(O)(O)=O")]

    phosphate_match = False
    for pattern in phosphate_patterns:
        if mol.HasSubstructMatch(pattern):
            phosphate_match = True
            break

    if not phosphate_match:
        return False, "No phosphate group found"

    # Look for sugar backbone patterns
    sugar_patterns = [Chem.MolFromSmarts("[OX2]r1[CX4]([OX2])[CX4]([OX2])[CX4]([OX2])1"),  # Ring
                      Chem.MolFromSmarts("[OX2][CX4]([OX2])[CX4]([OX2])[CX4]([OX2])"),  # Open chain
                      Chem.MolFromSmarts("[OX2]r1[CX4]([OX2])[CX4]([OX2])1")]  # 5-membered ring

    sugar_match = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(pattern):
            sugar_match = True
            break

    if not sugar_match:
        return False, "No sugar backbone found"

    # Check for phosphate-sugar connectivity
    phospho_sugar_pattern = Chem.MolFromSmarts("[P]~[OX2]~[CX4]")
    if not mol.HasSubstructMatch(phospho_sugar_pattern):
        return False, "Phosphate group not connected to sugar backbone"

    # Additional checks
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 3 or c_count > 20:
        return False, "Number of carbons outside the expected range"

    if o_count < 4 or o_count > 12:
        return False, "Number of oxygens outside the expected range"

    return True, "Contains a phosphate group attached to a sugar backbone"