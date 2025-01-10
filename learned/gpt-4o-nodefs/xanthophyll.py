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
        bool: True if the molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved pattern for long conjugated systems (flexible and extended)
    conjugated_pattern = Chem.MolFromSmarts("C=C(-C=C)*-C=C")  # Modify to better capture long conjugated chains
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_matches) < 3:
        return False, "No long conjugated carbon chain typical of xanthophylls found or chain too short"
    
    # Check for presence of various oxygen functionalities
    # Hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    # Carbonyl groups (=O)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    # Epoxides or cyclic ethers - more specific to capture rings
    epoxide_generic_pattern = Chem.MolFromSmarts("C1OC[C@@H]1")
    
    # Check for at least one of multiple oxygen functionalities
    oxygen_functionalities = any([
        mol.HasSubstructMatch(hydroxyl_pattern),
        mol.HasSubstructMatch(carbonyl_pattern),
        mol.HasSubstructMatch(epoxide_generic_pattern)
    ])
    
    if not oxygen_functionalities:
        return False, "No required oxygen functionalities found, necessary for xanthophylls"

    # Update molecular weight check to allow for variability in xanthophyll structures
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 400:  # Lower threshold to account for structure variability
        return False, f"Molecular weight {mol_weight} is lower than expected for typical xanthophyll"

    return True, "Contains long conjugated carbon chain and appropriate oxygen functionalities typical of xanthophylls"