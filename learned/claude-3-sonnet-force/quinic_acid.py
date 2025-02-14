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

    # Quinic acid core pattern with stereochemistry
    quinic_acid_pattern = Chem.MolFromSmarts("[OX2r6][C@@H]1[C@@](O)(C[C@H](O)[C@@H]([C@H]1O)C(O)=O)O")

    # Common substituents
    caffeoyl_pattern = Chem.MolFromSmarts("O=C/C=C/c:1:c:c:c:c:c:1O")
    feruloyl_pattern = Chem.MolFromSmarts("COc:1:c:c:c:c:c:1/C=C/C=O")
    coumaroyl_pattern = Chem.MolFromSmarts("O=C/C=C/c:1:c:c:c:c:c:1")

    # Check for quinic acid core
    if not mol.HasSubstructMatch(quinic_acid_pattern):
        return False, "No quinic acid core found"

    # Check for substituents
    sub_patterns = [caffeoyl_pattern, feruloyl_pattern, coumaroyl_pattern]
    sub_matches = [mol.GetSubstructMatches(p) for p in sub_patterns]
    num_subs = sum(len(matches) for matches in sub_matches)

    if num_subs == 0:
        return True, "Unsubstituted quinic acid"
    elif num_subs > 0:
        return True, f"Quinic acid with {num_subs} substituents"

    return False, "Not a quinic acid derivative"