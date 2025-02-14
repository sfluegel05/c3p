"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester as part of a carboxylic ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Acetate ester pattern: [-C(=O)OC(C)=O]
    # - This pattern ensures we are capturing a carbonyl single-bonded to an O, part of an acetic acid fragment
    acetate_ester_pattern = Chem.MolFromSmarts("C(=O)[O][CH3]")
    carboxylic_ester_pattern = Chem.MolFromSmarts("C(=O)O")  # Basic ester formation pattern

    if not mol.HasSubstructMatch(acetate_ester_pattern):
        return (False, "No acetate ester group found")

    if not mol.HasSubstructMatch(carboxylic_ester_pattern):
        return (False, "No ester functionality found")

    # Ensure we indeed focus on carboxylic esters with acetic acid component
    matches = mol.GetSubstructMatches(acetate_ester_pattern)
    for match in matches:
        if Chem.MolFromSmiles(smiles).GetAtomWithIdx(match[2]).GetSymbol() == 'C':
            # Ensure the methyl group is bonded properly to confirm it's an acetate ester
            return (True, "Contains correct acetate ester group")

    return (False, "Acetate ester group not part of acetic acid component")