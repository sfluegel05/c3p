"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are characterized by long conjugated carbon chains and the presence of specific oxygen functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for long conjugated chain pattern associated with carotenoids
    conjugated_pattern = Chem.MolFromSmarts("C=C(C)C=C(C)C=C(C)C=C(C)C=C(C)C=C(C)")  # Adjusted pattern for xanthophylls
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "No long conjugated carbon chain typical of xanthophylls found"
        
    # Check for presence of specific oxygen functionalities
    has_oxygen_functional_group = False
    # Check for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    # Check for carbonyl groups (=O)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    # Check for epoxide groups
    epoxide_pattern = Chem.MolFromSmarts("C1OC2(C1)C2")
    
    if mol.HasSubstructMatch(hydroxyl_pattern) or mol.HasSubstructMatch(carbonyl_pattern) or mol.HasSubstructMatch(epoxide_pattern):
        has_oxygen_functional_group = True
    
    if not has_oxygen_functional_group:
        return False, "No specific oxygen functionalities found, necessary for xanthophylls"

    # Check approximate molecular weight - xanthophylls are typically large molecules
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 500:
        return False, f"Molecular weight {mol_weight} too low for typical xanthophyll"

    return True, "Contains long conjugated carbon chain and specific oxygen functionalities typical of xanthophylls"