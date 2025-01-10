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
        bool: True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Attempt to define an improved core structure for dihydroagarofuran
    # This pattern improved only tentatively as research indicates these have complex polycyclic systems
    dihydroagarofuran_core_pattern = Chem.MolFromSmarts("O[C@@H]1[C@@H](OC)C[C@H]2[C@@H]1COC2")  # Improved hypothetical pattern

    if not mol.HasSubstructMatch(dihydroagarofuran_core_pattern):
        return False, "Core structure of dihydroagarofuran not found"

    # Check for ester groups - these should be more abundant in dihydroagarofuran
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 5:  # Based on the intricacy typical in the provided examples
        return False, f"Ester groups absent or fewer than expected, found {len(ester_matches)}"

    return True, "Molecule matches dihydroagarofuran sesquiterpenoid core structure with esters"