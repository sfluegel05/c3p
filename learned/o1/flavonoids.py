"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: flavonoids
"""

from rdkit import Chem

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a C6-C3-C6 skeleton consisting of two aromatic rings (A and B)
    connected via a heterocyclic ring (C) formed by a three-carbon bridge.
    This includes subclasses like flavones, isoflavones, neoflavonoids, chalcones, aurones,
    pterocarpans, flavonolignans, etc.

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

    # Core flavonoid scaffold: C6-C3-C6 with heterocyclic ring
    flavonoid_core = Chem.MolFromSmarts('c1cc(-c2ccc3occ(=O)c(-c4ccccc4)c3c2)ccc1')  # Flavone core

    # Isoflavone core: B ring at position 3
    isoflavone_core = Chem.MolFromSmarts('c1cc(-c2cc3oc(=O)cc(c3c2)c2ccccc2)ccc1')  # Isoflavone core

    # Flavanone core: Saturated heterocyclic ring
    flavanone_core = Chem.MolFromSmarts('c1cc(-c2ccc3occ(=O)cc3c2)ccc1')  # Flavanone core

    # Chalcone core: Open-chain flavonoid
    chalcone_core = Chem.MolFromSmarts('c1ccccc1C=CC(=O)c2ccccc2')  # Chalcone core

    # Dihydrochalcone core: Reduced chalcone
    dihydrochalcone_core = Chem.MolFromSmarts('c1ccccc1CCC(=O)c2ccccc2')  # Dihydrochalcone core

    # Aurone core
    aurone_core = Chem.MolFromSmarts('c1ccc(-c2oc(=O)c3ccccc3c2=O)cc1')  # Aurone core

    # Neoflavonoid core: B ring at position 4
    neoflavonoid_core = Chem.MolFromSmarts('c1cc(-c2ccc3oc(=O)cc(c3c2)c2ccccc2)ccc1')  # Neoflavonoid core

    # Pterocarpan core
    pterocarpan_core = Chem.MolFromSmarts('c1ccc2c(c1)c3c4ccc(O)cc4OC3CO2')  # Pterocarpan core

    # Flavonolignan core
    flavonolignan_core = Chem.MolFromSmarts('c1ccccc1C2OC3=CC=CC=C3OC2c4ccccc4')  # Flavonolignan core

    # Homoflavonoid core: C6-C3-C7 skeleton
    homoflavonoid_core = Chem.MolFromSmarts('c1ccccc1CC=CC(=O)c2ccccc2')  # Homoflavonoid core

    # Rotenoid core
    rotenoid_core = Chem.MolFromSmarts('c1cc2c(cc1)C(=O)c3c(c2)oc4ccccc34')  # Rotenoid core

    # Flavonoid oligomer (e.g., proanthocyanidins)
    flavonoid_oligomer_core = Chem.MolFromSmarts('c1cc(ccc1O)-c2c(O)cc(O)cc2O')  # Simplified pattern

    # List of patterns with their corresponding flavonoid subclass names
    patterns = [
        ('Flavonoid', flavonoid_core),
        ('Isoflavone', isoflavone_core),
        ('Flavanone', flavanone_core),
        ('Chalcone', chalcone_core),
        ('Dihydrochalcone', dihydrochalcone_core),
        ('Aurone', aurone_core),
        ('Neoflavonoid', neoflavonoid_core),
        ('Pterocarpan', pterocarpan_core),
        ('Flavonolignan', flavonolignan_core),
        ('Homoflavonoid', homoflavonoid_core),
        ('Rotenoid', rotenoid_core),
        ('Flavonoid Oligomer', flavonoid_oligomer_core),
    ]

    # Check for matches with flavonoid patterns
    for name, pattern in patterns:
        if pattern and mol.HasSubstructMatch(pattern):
            return True, f"Contains {name} core structure"

    # Additional check for C6-C3-C6 skeleton with heterocyclic ring
    # Count aromatic rings and heterocyclic rings
    ri = mol.GetRingInfo()
    num_aromatic_rings = 0
    num_hetero_rings = 0
    for ring in ri.BondRings():
        atoms_in_ring = [mol.GetBondWithIdx(idx).GetBeginAtomIdx() for idx in ring] + \
                        [mol.GetBondWithIdx(ring[-1]).GetEndAtomIdx()]
        is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in atoms_in_ring)
        contains_heteroatom = any(mol.GetAtomWithIdx(idx).GetAtomicNum() not in [6,1] for idx in atoms_in_ring)
        if is_aromatic:
            num_aromatic_rings += 1
        if contains_heteroatom:
            num_hetero_rings += 1

    if num_aromatic_rings >=2 and num_hetero_rings >=1:
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        # Flavonoids typically have 15 carbons (C15 skeleton)
        if carbon_count >= 15:
            return True, "Contains C6-C3-C6 skeleton with heterocyclic ring characteristic of flavonoids"

    return False, "No flavonoid core structure found"