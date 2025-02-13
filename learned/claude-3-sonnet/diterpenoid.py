"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: CHEBI:36655 diterpenoid

Diterpenoids are terpenoids derived from a diterpene, where the C20 skeleton of
the parent diterpene has been rearranged or modified by the removal of one or
more skeletal atoms (generally methyl groups).
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_diterpenoid(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight range (typical for diterpenoids)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.2f} outside typical diterpenoid range"

    # Check for diterpene backbone (20 carbon atoms)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Number of carbon atoms ({c_count}) not equal to 20"

    # Look for rearranged or modified skeletons
    sssr = Chem.GetSymmSSSR(mol)
    largest_ring = max(rings for rings in sssr)
    if len(largest_ring) < 4 or len(largest_ring) > 8:
        return False, "Diterpenoids typically have 4-8 membered ring systems"

    # Look for common diterpenoid substructures
    patterns = (
        "[C@H]1CC[C@@]2(C(=C)C)C[C@@H](O)[C@H](C)[C@]12C",  # labdane
        "[C@H]1C[C@@]2(C(=C)C)[C@@H](O)C[C@H](C)[C@]12C",   # clerodane
        "[C@H]1C[C@@H]2[C@@H](C)C[C@@H](C)[C@@]2(C)C1=C",   # pimarane
        "[C@H]1C[C@@H]2[C@@H](O)[C@H](C)C[C@@]2(C)C1=C",    # abietane
        "[C@H]1C[C@@H]2[C@@H](C)C[C@@H](C)[C@]2(C)C1=C",    # kaurane
    )

    for patt in patterns:
        patt_mol = Chem.MolFromSmarts(patt)
        if mol.HasSubstructMatch(patt_mol):
            return True, "Contains typical diterpenoid substructure"

    # No diterpenoid features found
    return False, "No diterpenoid features identified"