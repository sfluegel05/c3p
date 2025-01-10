"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a C6-C3-C6 carbon skeleton with a benzene ring
    attached to a heterocyclic benzopyran ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for flavonoid C6-C3-C6 structure with benzopyran ring (generic flavonoid pattern)
    # Basic flavonoid may be represented by a structure similar to 1-Benzopyran-4-one
    flavonoid_smart = "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]2:[#6](-[#6](=[#8])O2):[#6]:[#6]:1"  # Represents a benzopyranone (core of flavonoids)
    flavonoid_pattern = Chem.MolFromSmarts(flavonoid_smart)

    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "No flavonoid C6-C3-C6 backbone (benzopyran) found"

    return True, "Contains flavonoid C6-C3-C6 backbone with benzopyran structure"