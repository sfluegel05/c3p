"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: phosphatidic acid
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid is a derivative of glycerol where one hydroxyl group is esterified with phosphoric acid and the other two with fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is phosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define glycerol backbone pattern with ester linkages
    glycerol_pattern = Chem.MolFromSmarts("OCC(O[CX3](=O)O)CO[CX3](=O)O")
    # Define phosphoric acid ester pattern
    phosphoester_pattern = Chem.MolFromSmarts("COP(O)(O)=O")
    # Define fatty acid ester pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")

    # Match glycerol backbone including the phosphoric acid ester
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("OCC(*O)CO*")):
        return False, "No glycerol backbone pattern with ester linkages found"

    # Match phosphoric acid ester group
    if not mol.HasSubstructMatch(phosphoester_pattern):
        return False, "No phosphoric acid ester group found"

    # Match two esterified fatty acid groups
    ester_matches = list(mol.GetSubstructMatches(ester_pattern))
    if len(ester_matches) < 2:
        return False, "Less than 2 fatty acid ester groups found"

    return True, "Contains glycerol backbone with 2 fatty acid esters and a phosphoric acid ester"