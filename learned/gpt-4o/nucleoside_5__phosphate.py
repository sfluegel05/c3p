"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    Requires presence of a purine or pyrimidine base and a ribose or deoxyribose
    attached to a phosphate group at the 5' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify ribose or deoxyribose component
    ribose_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@@H](CO)O1")
    deoxyribose_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](CO)O1")

    has_ribose = mol.HasSubstructMatch(ribose_pattern)
    has_deoxyribose = mol.HasSubstructMatch(deoxyribose_pattern)
    
    if not has_ribose and not has_deoxyribose:
        return False, "No ribose or deoxyribose sugar found"

    # Identify purine or pyrimidine base
    purine_pattern = Chem.MolFromSmarts("n1cnc2n(cnc2c1)")
    pyrimidine_pattern = Chem.MolFromSmarts("c1cnc[nH]c1(=O)")

    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    
    if not (has_purine or has_pyrimidine):
        return False, "No purine or pyrimidine base found"

    # Identify phosphate group, ensuring it is attached to the sugar at the 5' position
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    
    if not phosphate_matches:
        return False, "No phosphate group found at C-5'"

    # Further checks can be added here to ensure the phosphate is monophosphate, diphosphate, etc.

    return True, "Contains a ribosyl or deoxyribosyl derivative of a purine/pyrimidine base with a 5'-phosphate"