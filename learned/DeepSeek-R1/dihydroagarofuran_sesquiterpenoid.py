"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES.
    These compounds feature a bridged bicyclic skeleton characteristic of dihydroagarofuran derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if matches criteria, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Revised core pattern: bridged bicyclic system without stereochemistry or specific substituents
    # This SMARTS aims to capture the bicyclo[5.3.0]decane-like skeleton
    core_pattern = Chem.MolFromSmarts("C12CC3CC(C1)C2C3")
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing dihydroagarofuran core structure"
    
    # Additional check for sesquiterpenoid: 15 carbons in the parent skeleton
    # Note: This is a simplified check and may not be accurate for all derivatives
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Insufficient carbons for sesquiterpenoid (minimum 15)"
    
    return True, "Contains dihydroagarofuran skeleton with sesquiterpenoid characteristics"