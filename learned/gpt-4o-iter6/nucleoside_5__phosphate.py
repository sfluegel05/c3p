"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate consists of a ribose or deoxyribose sugar linked to a purine or pyrimidine base,
    with phosphorylation at the C-5 position of the sugar.

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

    # Define sugar pattern for both ribose and deoxyribose
    sugar_pattern = Chem.MolFromSmarts("C1([C@@H]2C(O)[C@@H](O)CO2)O1")  # General 5-membered sugar ring

    # Generalized phosphate group pattern for mono-, di-, tri-, tetra-phosphate
    phosphate_pattern = Chem.MolFromSmarts("COP(O)(=O)(O)O")  # Attaching phosphate to ribose

    # Simplified purine and pyrimidine patterns
    purine_pattern = Chem.MolFromSmarts("c1ncnc2ncnc12")
    pyrimidine_pattern = Chem.MolFromSmarts("c1[nH]cnc(N)c1=O")

    # Check for sugar pattern matches in the molecule
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "Does not contain a ribose or deoxyribose sugar structure"

    # Check for phosphate group patterns
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found attached at 5' position of the sugar"

    # Check for either purine or pyrimidine base
    base_matches = (
        mol.HasSubstructMatch(purine_pattern) or 
        mol.HasSubstructMatch(pyrimidine_pattern)
    )
    if not base_matches:
        return False, "No recognized purine or pyrimidine base structure found"

    return True, "Contains a ribosyl or deoxyribosyl sugar, phosphate at 5', and a nucleobase"