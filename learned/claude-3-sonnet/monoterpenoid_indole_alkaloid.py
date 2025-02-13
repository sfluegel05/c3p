"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: CHEBI:24049 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for indole ring system
    indole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)nc(n2)C")
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole ring system found"

    # Check for terpene unit (fused ring system with 10 carbon atoms)
    terpene_pattern = Chem.MolFromSmarts("C1CCCCCCCCC1")
    if not mol.HasSubstructMatch(terpene_pattern):
        return False, "No terpene unit found"

    # Check for nitrogen atom and heavy atom count
    if not any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "No nitrogen atom found"

    heavy_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1)
    if heavy_atom_count < 20 or heavy_atom_count > 40:
        return False, "Heavy atom count outside the expected range"

    return True, "Contains an indole ring system and a terpene unit, with a nitrogen atom and appropriate heavy atom count"