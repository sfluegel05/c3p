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

    # Define SMARTS patterns for core flavonoid structures

    # Flavone core (2-phenylchromen-4-one)
    flavone_pattern = Chem.MolFromSmarts('O=C1C=CC2=CC=CC=C2Oc3c1cccc3')

    # Isoflavone core (3-phenylchromen-4-one)
    isoflavone_pattern = Chem.MolFromSmarts('O=C1C=CC2=CC=CC=C2C=Cc3c1cccc3')

    # Flavanone core (2-phenylchroman-4-one)
    flavanone_pattern = Chem.MolFromSmarts('O=C1CCC2=CC=CC=C2Oc3c1cccc3')

    # Chalcone core (1,3-diphenylprop-2-en-1-one)
    chalcone_pattern = Chem.MolFromSmarts('O=CC=CCc1ccccc1c2ccccc2')

    # Dihydrochalcone core (1,3-diphenylpropan-1-one)
    dihydrochalcone_pattern = Chem.MolFromSmarts('O=CCCc1ccccc1c2ccccc2')

    # Aurone core (2-benzylidenebenzofuran-3-one)
    aurone_pattern = Chem.MolFromSmarts('O=C1OC=CC1=CC=c2ccccc2')

    # Flavanol core (2-phenyl-3,4-dihydro-2H-chromen-3-ol)
    flavanol_pattern = Chem.MolFromSmarts('Oc1cc(O)ccc1C2OCc3ccccc23')

    # Anthocyanidin core (2-phenylbenzopyrylium)
    anthocyanidin_pattern = Chem.MolFromSmarts('[O+]1c2ccccc2Oc3ccccc13')

    # Pterocarpan core (benzopyran fused with dihydrofuran ring)
    pterocarpan_pattern = Chem.MolFromSmarts('O1CC2Oc3ccccc3C4=C2C=CC=C14')

    # List of patterns to check
    flavonoid_patterns = [
        (flavone_pattern, "Flavone core detected"),
        (isoflavone_pattern, "Isoflavone core detected"),
        (flavanone_pattern, "Flavanone core detected"),
        (chalcone_pattern, "Chalcone core detected"),
        (dihydrochalcone_pattern, "Dihydrochalcone core detected"),
        (aurone_pattern, "Aurone core detected"),
        (flavanol_pattern, "Flavanol core detected"),
        (anthocyanidin_pattern, "Anthocyanidin core detected"),
        (pterocarpan_pattern, "Pterocarpan core detected"),
    ]

    # Check for each pattern
    for pattern, reason in flavonoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, reason

    # Check for C6-C3-C6 skeleton (general flavonoid structure)
    c6c3c6_pattern = Chem.MolFromSmarts('c1ccc(cc1)CCCc2ccccc2')
    if mol.HasSubstructMatch(c6c3c6_pattern):
        return True, "General C6-C3-C6 flavonoid skeleton detected"

    return False, "No flavonoid core structure detected"