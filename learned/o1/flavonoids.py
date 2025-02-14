"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: flavonoids
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a C6-C3-C6 skeleton consisting of two aromatic rings (A and B)
    connected via a three-carbon bridge that may form a heterocyclic ring (C).
    
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

    # Define SMARTS patterns for flavonoid subclasses

    # Flavone core: 2-phenylchromen-4-one
    flavone_pattern = Chem.MolFromSmarts('O=C1C=CC2=CC=CC=C2O1')  # Flavone core
    # Isoflavone core: 3-phenylchromen-4-one
    isoflavone_pattern = Chem.MolFromSmarts('O=C1C=CC(Oc2ccccc2)=C1')  # Isoflavone core
    # Flavanone core: 2,3-dihydroflavone
    flavanone_pattern = Chem.MolFromSmarts('O=C1CCc2ccccc2O1')  # Flavanone core
    # Chalcone core: 1,3-diphenylprop-2-en-1-one
    chalcone_pattern = Chem.MolFromSmarts('O=CC=CCc1ccccc1')  # Chalcone core
    # Aurone core: 2-benzylidenebenzofuran-3(2H)-one
    aurone_pattern = Chem.MolFromSmarts('O=C1C=CC(=O)c2ccccc12')  # Aurone core
    # Flavanol core: 3-hydroxyflavanone
    flavanol_pattern = Chem.MolFromSmarts('O1CCC(CC1)C2=CC=CC=C2')  # Flavanol core

    patterns = [
        ('Flavone', flavone_pattern),
        ('Isoflavone', isoflavone_pattern),
        ('Flavanone', flavanone_pattern),
        ('Chalcone', chalcone_pattern),
        ('Aurone', aurone_pattern),
        ('Flavanol', flavanol_pattern)
    ]
    
    for name, pattern in patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return True, f"Contains {name} core structure"

    # Check for general C6-C3-C6 skeleton
    # Two aromatic rings connected via a three-carbon chain or heterocycle
    c6_c3_c6_pattern = Chem.MolFromSmarts('c1ccc(-c2ccc(-c3ccccc3)cc2)cc1')  # General C6-C3-C6 pattern
    if mol.HasSubstructMatch(c6_c3_c6_pattern):
        return True, "Contains C6-C3-C6 skeleton characteristic of flavonoids"

    return False, "No flavonoid core structure found"