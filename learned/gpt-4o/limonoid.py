"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids potentially containing a 4,4,8-trimethyl-17-furanylsteroid skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for furan ring; limonoids often have one but this is not exclusive
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    furan_match = mol.HasSubstructMatch(furan_pattern)

    # Count oxygen atoms; limonoids are highly oxygenated
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, f"Insufficient oxygenation, found {o_count} oxygens"

    # Count carbon atoms; consistent with triterpenoids range
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 40:
        return False, f"Unusual triterpenoid carbon count, found {c_count} carbons"

    # Check if it is a polycyclic compound
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 4:
        return False, f"Too few rings, found {n_rings}"

    # Even if a furan isn't present, ensure sufficient heteroatoms and ring complexity
    if furan_match:
        return True, "Contains a furan ring and matches typical structural features of limonoids"
    else:
        return True, "No furan ring, but matches other limonoid structural features (high oxygenation, polycyclic skeleton)"

# This function aims to classify a SMILES string as representing a limonoid by considering
# several criteria common in limonoid structures but offers flexibility to capture structural diversity.