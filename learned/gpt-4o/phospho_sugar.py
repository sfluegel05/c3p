"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is a monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

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

    # Generic sugar pattern (5 or 6-membered ring with oxygens)
    sugar_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](CO)[C@@H](O)[C@@H]1O | OC[C@H]1O[C@H](O)[C@H](O)[C@@H]1O | OC1C(O)C(O)C(O)CO1 | OC1C(O)C(O)C(O)C1")
    if not sugar_pattern or not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar ring structure found"
    
    # Phosphate ester linkage: -O-P(=O)(O)O
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    if not phosphate_pattern or not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate ester linkage found"

    return True, "Contains sugar ring structure with a phosphate ester linkage"