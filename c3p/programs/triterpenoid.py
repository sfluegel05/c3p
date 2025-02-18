"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: CHEBI:26870 triterpenoid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    A triterpenoid is derived from a triterpene, typically containing 30 carbons,
    but may have modifications such as rearrangements or loss of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 24 or c_count > 33:
        return False, f"Carbon count {c_count} not in typical triterpenoid range (24-33 carbons)"
    
    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, likely not a triterpenoid"

    # Calculate number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 4:
        return False, f"Only {num_rings} rings found, triterpenoids typically have multiple ring structures"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight {mol_wt:.2f} Da too low for triterpenoid"

    return True, "Molecule meets criteria for a triterpenoid (carbon count, oxygen atoms, ring structures)"