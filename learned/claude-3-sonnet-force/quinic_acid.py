"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: CHEBI:18198 quinic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid based on its SMILES string.
    A quinic acid is a cyclitol carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for quinic acid core pattern
    quinic_acid_core = Chem.MolFromSmarts("[OX2r6][CX4r6][CX4r6][CX4r6][CX4r6][CX4r6]1[OX2H1,OX1H0-]")
    if not mol.HasSubstructMatch(quinic_acid_core):
        return False, "No quinic acid core found"

    # Look for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Look for common substituents (caffeoyl, feruloyl, coumaroyl)
    caffeoyl_pattern = Chem.MolFromSmarts("O=C/C=C/c:1:c:c:c:c:c:1O")
    feruloyl_pattern = Chem.MolFromSmarts("COc:1:c:c:c:c:c:1/C=C/C=O")
    coumaroyl_pattern = Chem.MolFromSmarts("O=C/C=C/c:1:c:c:c:c:c:1")

    sub_patterns = [caffeoyl_pattern, feruloyl_pattern, coumaroyl_pattern]
    sub_matches = [mol.GetSubstructMatches(p) for p in sub_patterns]
    num_subs = sum(len(matches) for matches in sub_matches)

    if num_subs > 3:
        return False, "Too many substituents for quinic acid"

    return True, "Contains quinic acid core with carboxylic acid group"