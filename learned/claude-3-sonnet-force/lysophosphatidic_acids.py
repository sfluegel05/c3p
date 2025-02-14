"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: CHEBI:26676 lysophosphatidic acids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid is a monoacylglycerol phosphate, with one or more fatty acid chains attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol moiety (C-C-C with 2 or 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X3,CH2X4][CHX3,CHX4][CH2X3,CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol moiety found"
        
    # Look for phosphate group (-O-P(=O)(O)O)
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for fatty acid chains (long carbon chains)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "No significant fatty acid chains found"

    # Check molecular weight - lysophosphatidic acids typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for lysophosphatidic acid"

    return True, "Contains glycerol moiety with one or more fatty acid chains and a phosphate group"