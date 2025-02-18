"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: CHEBI:23899 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    Diterpenoids are derived from a C20 diterpene skeleton, which may be modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check carbon count (allowing for some modification)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 22):
        return False, f"Carbon count {c_count} outside diterpenoid range (18-22)"

    # Check if it's a terpenoid (at least 4 isoprene units)
    terpene_units = rdMolDescriptors.CalcNumTerpenes(mol)
    if terpene_units < 4:
        return False, f"Only {terpene_units} terpene units, need ≥4 for diterpenoid"

    return True, f"{c_count} carbons and {terpene_units} terpene units indicate diterpenoid"