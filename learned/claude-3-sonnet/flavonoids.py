"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: CHEBI:25107 flavonoids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a C15 or C16 skeleton with a phenyl-substituted 1-phenylpropane structure,
    or a condensed C6-C3 lignan precursor.

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

    # Define flavonoid core structure patterns
    flavonoid_core_patterns = [
        Chem.MolFromSmarts("[c1c(cc2c(c1)OCc3c(c2)cccc3)O]"),   # Flavone
        Chem.MolFromSmarts("[c1c(cc2c(c1)OCC3=CC(=O)Oc4c3cccc4)O]"),  # Flavonol
        Chem.MolFromSmarts("[c1c(cc2c(c1)OCC3=CC(=O)c4c(O)cccc4O3)O]"),  # Flavanone
        Chem.MolFromSmarts("[c1c(c2c(cc1O)OCC3=CC(=O)c4c(O)cccc4O3)O]"),  # Flavan-3-ol
        Chem.MolFromSmarts("[c1c(c2c(cc1O)OCC3=C(O)c4c(cccc4O3)O)O]"),  # Flavan-4-ol
        Chem.MolFromSmarts("[c1c(c2c(cc1O)OCC3=C(O)c4c(cc(O)cc4O3)O)O]"),  # Flavan-3,4-diol
        Chem.MolFromSmarts("[c1c(c2c(cc1O)OCC3=CC(=O)c4c(ccc(O)c4O3)O)O]"),  # Isoflavone
        Chem.MolFromSmarts("[c1c(c2c(cc1O)Oc3c(c2)cccc3)O]"),  # Chalcone
        Chem.MolFromSmarts("[c1c(c2c(cc1O)Oc3c(c2)c(O)ccc3)O]"),  # Aurone
        Chem.MolFromSmarts("[c1c(c2c(cc1O)Oc3c(c2)cccc3O)O]"),  # Dihydrochalcone
        Chem.MolFromSmarts("[c1c(c2c(cc1O)Oc3c(c2)cccc3)CO]"),  # Neoflavonoid
    ]

    # Check for flavonoid core structure
    for pattern in flavonoid_core_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule contains a flavonoid core structure"

    # Check for lignan precursor
    lignan_pattern = Chem.MolFromSmarts("[c1c(ccc2c1Oc3ccccc3O2)O]")
    if mol.HasSubstructMatch(lignan_pattern):
        return True, "Molecule contains a condensed C6-C3 lignan precursor"

    # Check for additional structural rules
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    num_aromatic_rings = sum(1 for ring in ring_info.AtomRings() if ring.IsAromatic())

    # Flavonoids should have at least one aromatic ring
    if num_aromatic_rings == 0:
        return False, "No aromatic rings found"

    # Flavonoids should have a total of 2-4 rings
    if num_rings < 2 or num_rings > 4:
        return False, f"Found {num_rings} rings, but flavonoids should have 2-4 rings"

    # Flavonoids should have at least one oxygen atom
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygens == 0:
        return False, "No oxygen atoms found"

    # Flavonoids should have a molecular weight between 200-800 Da
    mol_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 200 or mol_weight > 800:
        return False, "Molecular weight outside the typical range for flavonoids"

    return True, "Molecule satisfies the structural rules for flavonoids"