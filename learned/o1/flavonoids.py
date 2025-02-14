"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: flavonoids
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a C6-C3-C6 skeleton consisting of two aromatic rings (A and B)
    connected via a heterocyclic ring (C) formed by a three-carbon bridge.
    This includes subclasses like flavones, isoflavones, neoflavonoids, chalcones, aurones, etc.

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

    # Define SMARTS patterns for different flavonoid subclasses

    # Flavone core: 2-phenylchromen-4-one
    flavone_pattern = Chem.MolFromSmarts('O=C1C=CC=CC1=O-c2ccccc2')

    # Isoflavone core: 3-phenylchromen-4-one
    isoflavone_pattern = Chem.MolFromSmarts('O=C1C=CC=CC1=COc2ccccc2')

    # Flavanone core: 2,3-dihydroflavone
    flavanone_pattern = Chem.MolFromSmarts('O=C1CC=CC=CC1=O-c2ccccc2')

    # Flavanol core: 3-hydroxyflavanone
    flavanol_pattern = Chem.MolFromSmarts('O=C1CC(OC2=CC=CC=C2)C=CC1=O')

    # Chalcone core: 1,3-diphenylprop-2-en-1-one
    chalcone_pattern = Chem.MolFromSmarts('O=CC=CCc1ccccc1')

    # Aurone core: 2-benzylidenebenzofuran-3(2H)-one
    aurone_pattern = Chem.MolFromSmarts('O=C1OC=CC2=CC=CC=C12')

    # Dihydrochalcone core
    dihydrochalcone_pattern = Chem.MolFromSmarts('O=CCCc1ccccc1')

    # Neoflavonoid core
    neoflavonoid_pattern = Chem.MolFromSmarts('c1cc(O)ccc1C2=COC(=O)c3ccccc23')

    # Pterocarpan core
    pterocarpan_pattern = Chem.MolFromSmarts('c1ccc2c(c1)C3Oc4ccccc4C(C2O3)')

    # Flavonolignan core (e.g., silybin structure)
    flavonolignan_pattern = Chem.MolFromSmarts('C1COc2cc3c(cc2O1)-c1ccccc1O3')

    # List of patterns with their corresponding flavonoid subclass names
    patterns = [
        ('Flavone', flavone_pattern),
        ('Isoflavone', isoflavone_pattern),
        ('Flavanone', flavanone_pattern),
        ('Flavanol', flavanol_pattern),
        ('Chalcone', chalcone_pattern),
        ('Aurone', aurone_pattern),
        ('Dihydrochalcone', dihydrochalcone_pattern),
        ('Neoflavonoid', neoflavonoid_pattern),
        ('Pterocarpan', pterocarpan_pattern),
        ('Flavonolignan', flavonolignan_pattern),
    ]

    # Check for matches with flavonoid patterns
    for name, pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {name} core structure"

    # Additional check for C6-C3-C6 skeleton with heterocyclic ring
    # Identify rings in the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    aromatic_rings = [ring for ring in atom_rings if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]

    if len(aromatic_rings) >= 2:
        # Check for a three-carbon bridge connecting two aromatic rings
        for bond in mol.GetBonds():
            if bond.GetBondTypeAsDouble() == 1:
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                if begin_atom.GetIsAromatic() and end_atom.GetIsAromatic():
                    path = Chem.GetShortestPath(mol, begin_atom.GetIdx(), end_atom.GetIdx())
                    # Exclude direct connections between aromatic atoms (should be part of the same ring)
                    if len(path) == 4:
                        # Check if the path includes three connecting atoms (C3 bridge)
                        bridge_atoms = path[1:-1]
                        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in bridge_atoms):
                            return True, "Contains C6-C3-C6 skeleton characteristic of flavonoids"

    return False, "No flavonoid core structure found"