"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are characterized by long conjugated carbon chains and various oxygen-containing functional groups.

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

    # General pattern for long conjugated chains (flexible pattern)
    conjugated_pattern = Chem.MolFromSmarts("C=C(-C)*C=C")  # Extended pattern for conjugated systems
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "No long conjugated carbon chain typical of xanthophylls found"

    # Check for presence of various oxygen functionalities
    # Hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    # Carbonyl groups (=O)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    # Cyclic ethers, epoxides, or other oxygen-containing rings
    epoxide_generic_pattern = Chem.MolFromSmarts("C1OC=C1")
    
    # Check for multiple oxygen functionalities
    oxygen_functionalities = any([
        mol.HasSubstructMatch(hydroxyl_pattern),
        mol.HasSubstructMatch(carbonyl_pattern),
        mol.HasSubstructMatch(epoxide_generic_pattern)
    ])
    
    if not oxygen_functionalities:
        return False, "No required oxygen functionalities found, necessary for xanthophylls"

    # Check molecular weight appropriate for large, complex xanthophyll structures
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 500:
        return False, f"Molecular weight {mol_weight} is lower than expected for typical xanthophyll"

    return True, "Contains long conjugated carbon chain and appropriate oxygen functionalities typical of xanthophylls"