"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid or its derivative based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic quinic acid pattern: a cyclohexane with multiple hydroxy groups and one carboxylic acid
    # Allow flexibility in stereochemistry with '[C@H?' to account for stereocenters variability
    quinic_acid_base_pattern = Chem.MolFromSmarts("OC1[C@H?](O)[C@H?](O)[C@H?](O)[C@H?](O)C1C(O)=O")
    
    # Check if molecule contains the quinic acid core
    if not mol.HasSubstructMatch(quinic_acid_base_pattern):
        return False, "Mismatch in quinic acid core structure"

    # Additional checks for common quinic acid derivatives
    # Check for presence of common ester linkages such as caffeoyl or other similar substitutions
    ester_pattern = Chem.MolFromSmarts("O=C(O)C=Cc1ccc(O)c(O)c1")
    if mol.HasSubstructMatch(ester_pattern):
        return True, "Matches quinic acid core structure with common ester substitutions"

    return True, "Matches quinic acid core structure"