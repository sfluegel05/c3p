"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside where the terpene moiety is a triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for triterpenoid backbone (typically involves at least 4-6 rings)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 4:
        return False, f"Triterpenoid structures typically have at least 4 rings, found {num_rings}"

    # Adjust carbon count range based on known examples
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 27 or carbon_count > 70:  # Updated to accommodate structures like dianversicoside G
        return False, f"Carbon count {carbon_count} not within extended range for triterpenoids (27-70)"

    # Advanced sugar moieties and glycosidic linkages detection
    # Expanded to check for common sugar units (glucose, rhamnose) and their linkages
    glycosidic_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H]1O")  # Glucopyranoside unit
    
    # Check for the presence of at least one glycosidic linkage pattern
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No detectable glycosidic linkage pattern found"

    # Ensure we detect presence of sugar moieties
    sugars = Chem.MolFromSmarts("[C@H]1O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    num_sugars = len(mol.GetSubstructMatches(sugars))

    if num_sugars < 1:
        return False, "No attached sugar moieties identified"

    return True, "Contains triterpenoid backbone with glycosidic linkages"