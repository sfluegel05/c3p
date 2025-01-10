"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for nucleobase pattern (purine or pyrimidine)
    purine_pattern = Chem.MolFromSmarts('n1cnc2c1ncnc2N')
    pyrimidine_pattern = Chem.MolFromSmarts('n1ccncn1')
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No nucleobase (purine or pyrimidine) found"
    
    # Look for sugar pattern (ribose or deoxyribose)
    ribose_pattern = Chem.MolFromSmarts('C1C(O)C(O)C(O)C1O')
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose or deoxyribose sugar found"

    # Look for phosphate group(s) (-O-P(=O)(O)-)
    phosphate_pattern = Chem.MolFromSmarts('O=P(O)(O)O')
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate group(s) found"

    return True, "Molecule contains a nucleoside structure with one or more phosphate groups"