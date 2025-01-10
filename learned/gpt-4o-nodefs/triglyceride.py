"""
Classifies: CHEBI:17855 triglyceride
"""
from rdkit import Chem

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride is defined as a glycerol esterified with three fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triglyceride, False otherwise
        str: Reason for classification
    """
    
    # 1. Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Look for esterified glycerol pattern 
    # (Chiral Centers may complicate pattern detection)
    glycerol_ester_pattern = Chem.MolFromSmarts("COC(=O)[CH](COC(=O))COC(=O)")
    if not mol.HasSubstructMatch(glycerol_ester_pattern):
        return False, "No glycerol ester backbone found"

    # 3. Check for three ester linkages
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 3:
        return False, f"Found {len(ester_matches)} ester groups, need at least 3"
    
    # 4. Check for sufficient fatty acid chain length through rotatable bonds
    n_rotatable = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains not long enough to represent typical fatty acids"
    
    return True, "Contains glycerol ester backbone with at least three ester linkages and fatty acid chains"