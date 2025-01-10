"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: Flavonoids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are polyphenolic compounds consisting of two phenyl rings (A and B)
    connected by a three-carbon bridge that may form a third heterocyclic ring (C),
    resulting in a C6-C3-C6 skeleton.

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

    # Define more general SMARTS patterns for core flavonoid structures

    # General flavonoid core (C6-C3-C6 skeleton)
    flavonoid_core_pattern = Chem.MolFromSmarts('[cR6]1[c,c]ccc[c,c]1[$([CH2][CH2][CH2]),$([CH]=[CH][CH]=),$([CH2][CH]=[CH2]),$([CH]=[CH][CH2])][cR6]2[c,c]ccc[c,c]2')

    # Flavone and flavonol core (2-phenylchromen-4-one and derivatives)
    flavone_pattern = Chem.MolFromSmarts('O=C1C=CC(=O)c2c(O1)cccc2')

    # Isoflavone core (3-phenylchromen-4-one)
    isoflavone_pattern = Chem.MolFromSmarts('O=C1C=CC(=O)c2c(C1=CC=CC2)cccc2')

    # Flavanone core (2-phenylchroman-4-one)
    flavanone_pattern = Chem.MolFromSmarts('O=C1CCc2c(O1)cccc2')

    # Chalcone core (1,3-diphenylprop-2-en-1-one)
    chalcone_pattern = Chem.MolFromSmarts('O=CC=CCc1ccccc1')

    # Dihydrochalcone core (1,3-diphenylpropan-1-one)
    dihydrochalcone_pattern = Chem.MolFromSmarts('O=CCCc1ccccc1')

    # Aurone core (2-benzylidenebenzofuran-3-one)
    aurone_pattern = Chem.MolFromSmarts('O=C1OC=CC1=CCc2ccccc2')

    # Flavanol core (2-phenyl-3,4-dihydro-2H-chromen-3-ol)
    flavanol_pattern = Chem.MolFromSmarts('O[C@@H]1CCc2ccccc2O1')

    # Anthocyanidin core (2-phenylbenzopyrylium)
    anthocyanidin_pattern = Chem.MolFromSmarts('[O+]1c2ccccc2Oc3ccccc13')

    # Pterocarpan core (benzopyran fused with dihydrofuran ring)
    pterocarpan_pattern = Chem.MolFromSmarts('O1CC2Oc3ccccc3C4=C2C=CC=C14')

    # Neoflavonoid core (4-arylchromene)
    neoflavonoid_pattern = Chem.MolFromSmarts('c1ccc2Oc3ccccc3C=Cc2c1')

    # Include flavonolignan core
    flavonolignan_pattern = Chem.MolFromSmarts('c1ccccc1C2CC(OC3=CC=CC=C3O2)c4ccccc4')

    # List of patterns to check
    flavonoid_patterns = [
        (flavonoid_core_pattern, "General C6-C3-C6 flavonoid skeleton detected"),
        (flavone_pattern, "Flavone or flavonol core detected"),
        (isoflavone_pattern, "Isoflavone core detected"),
        (flavanone_pattern, "Flavanone core detected"),
        (chalcone_pattern, "Chalcone core detected"),
        (dihydrochalcone_pattern, "Dihydrochalcone core detected"),
        (aurone_pattern, "Aurone core detected"),
        (flavanol_pattern, "Flavanol core detected"),
        (anthocyanidin_pattern, "Anthocyanidin core detected"),
        (pterocarpan_pattern, "Pterocarpan core detected"),
        (neoflavonoid_pattern, "Neoflavonoid core detected"),
        (flavonolignan_pattern, "Flavonolignan core detected"),
    ]

    # Check for each pattern
    for pattern, reason in flavonoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, reason

    return False, "No flavonoid core structure detected"