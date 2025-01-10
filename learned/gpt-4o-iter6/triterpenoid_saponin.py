"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is defined as a terpene glycoside where the terpene moiety is a triterpenoid.

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
    
    # Look for triterpenoid backbone pattern
    # This is a complex pattern; approximate through number of rings and carbon count
    if rdMolDescriptors.CalcNumRings(mol) < 4:
        return False, "Triterpenoid structures commonly have at least 4 rings"

    # Further check if backbone resembles typical triterpene structure (30 carbon atoms, potentially more depending on the backbone)
    # Simplified characteristic of triterpenoid
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 27 or carbon_count > 36:  # Triterpenoid core typically has 30 carbons
        return False, f"Carbon count {carbon_count} not within typical range for triterpenoids"

    # Look for glycosidic linkage patterns
    # Using a generic pattern for sugar connections (O-C glycosidic bond)
    glycosidic_pattern = Chem.MolFromSmarts("C-O-C")
    
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage detected"

    # Count the number of potential sugar moieties
    sugars = Chem.MolFromSmarts("[C@H]1O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]1O |1:2,3,4,5,6|")
    num_sugars = mol.GetSubstructMatches(sugars)

    if len(num_sugars) == 0:
        return False, "No attached sugar moieties found"

    return True, "Contains triterpenoid backbone with glycosidic linkages"