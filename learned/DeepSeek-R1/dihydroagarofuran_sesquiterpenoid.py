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
    
    # Define the core dihydroagarofuran skeleton pattern
    # This SMARTS captures the bridged bicyclic system with key oxygen placements
    core_pattern = Chem.MolFromSmarts("[C@]12[C@@H](O)[C@H](C)[C@@H](C)C[C@H]1C[C@H]2C")
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing dihydroagarofuran core structure"
    
    return True, "Contains dihydroagarofuran skeleton with sesquiterpenoid characteristics"