"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: Flavonoids
"""
from rdkit import Chem

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are polyphenolic compounds consisting of two phenyl rings (A and B)
    connected by a three-carbon bridge that may form a third heterocyclic ring (C),
    resulting in a C6-C3-C6 skeleton.

    This function checks for the presence of the flavonoid core structure,
    allowing for various substitutions and modifications.

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

    # General flavonoid core pattern (C6-C3-C6 skeleton)
    # This pattern represents two aromatic rings (A and B) connected via
    # a heterocyclic ring (ring C), which may have various substituents.

    # Flavone/Flavonol core (allows substitutions at any position)
    flavone_core = Chem.MolFromSmarts('O=C1C=CC(=C[C;R]2)Oc3ccccc3C1c4ccccc24')

    # Flavanone core
    flavanone_core = Chem.MolFromSmarts('O=C1CC[C;R]2Oc3ccccc3C1c4ccccc24')

    # Isoflavone core
    isoflavone_core = Chem.MolFromSmarts('O=C1C=CC(=C[C;R]2)Oc3ccccc3C1c4ccccc24')

    # Chalcone core (open chain)
    chalcone_core = Chem.MolFromSmarts('O=CC=CCc1ccccc1')

    # Dihydrochalcone core
    dihydrochalcone_core = Chem.MolFromSmarts('O=CCC[C;R]1c2ccccc2cc1')

    # Aurone core
    aurone_core = Chem.MolFromSmarts('O=C1Oc2ccccc2C=C1c3ccccc3')

    # Flavanol core (with hydroxyl group at position 3)
    flavanol_core = Chem.MolFromSmarts('O[C@H]1C[C;R]2Oc3ccccc3C1c4ccccc24')

    # Anthocyanidin core (flavylium cation)
    anthocyanidin_core = Chem.MolFromSmarts('[O+]C1=C([C;R]2)Oc3ccccc3C=C1c4ccccc24')

    # Neoflavonoid core
    neoflavonoid_core = Chem.MolFromSmarts('c1ccc2Oc3ccccc3C=C2c1')

    # Pterocarpan core
    pterocarpan_core = Chem.MolFromSmarts('O1CC2Oc3ccccc3C4=C2C=CC=C14')

    # Flavonolignan core
    flavonolignan_core = Chem.MolFromSmarts('c1ccccc1[C;R]2CC(OC3=CC=CC=C3O2)c4ccccc4')

    # List of patterns to check with their corresponding reasons
    flavonoid_patterns = [
        (flavone_core, "Flavone or flavonol core detected"),
        (flavanone_core, "Flavanone core detected"),
        (isoflavone_core, "Isoflavone core detected"),
        (chalcone_core, "Chalcone core detected"),
        (dihydrochalcone_core, "Dihydrochalcone core detected"),
        (aurone_core, "Aurone core detected"),
        (flavanol_core, "Flavanol core detected"),
        (anthocyanidin_core, "Anthocyanidin core detected"),
        (neoflavonoid_core, "Neoflavonoid core detected"),
        (pterocarpan_core, "Pterocarpan core detected"),
        (flavonolignan_core, "Flavonolignan core detected"),
    ]

    # Check for each pattern
    for pattern, reason in flavonoid_patterns:
        if pattern is None:
            continue  # Skip invalid pattern
        if mol.HasSubstructMatch(pattern):
            return True, reason

    # As an additional check, look for the general C6-C3-C6 skeleton
    # This pattern matches two aromatic rings connected by three atoms (which may be part of a ring)

    general_core = Chem.MolFromSmarts('[a]1aaaaa1-[C;!R](-[C;!R]-[C;!R])-a2aaaaa2')
    if mol.HasSubstructMatch(general_core):
        return True, "General flavonoid core (C6-C3-C6 skeleton) detected"

    return False, "No flavonoid core structure detected"