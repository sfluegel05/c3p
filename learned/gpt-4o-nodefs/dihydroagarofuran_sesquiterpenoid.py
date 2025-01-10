"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    The classification is based on known core structures and functional group features typical of this class.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core structure pattern for dihydroagarofuran, adjusting for complexity and stereochemistry
    # Tentative pattern based on analysis and correction from previous version
    dihydroagarofuran_core_pattern = Chem.MolFromSmarts("OC1[CH]2OC(=O)[C@@]3(C)C[C@H](OC(=O)C)C[C@@H]3C[C@H]12")  

    if not mol.HasSubstructMatch(dihydroagarofuran_core_pattern):
        return False, "Core structure of dihydroagarofuran not found"

    # Check for multiple ester groups as indicators of high esterification
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 5:  # Based on the high degree of esterification typical in the provided examples
        return False, f"Ester groups absent or fewer than expected, found {len(ester_matches)}"

    return True, "Molecule matches dihydroagarofuran sesquiterpenoid core structure with esters"