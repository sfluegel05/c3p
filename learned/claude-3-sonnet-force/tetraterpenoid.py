"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:35625 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid is a terpenoid derived from a tetraterpene (C40 skeleton),
    where the parent tetraterpene skeleton may have been rearranged or modified
    by the removal of one or more skeletal atoms (generally methyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for tetraterpene skeleton (C40)
    num_atoms = mol.GetNumAtoms()
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 40:
        return False, f"Does not have a C40 skeleton (found {num_carbons} carbons)"

    # Check for terpenoid-like structure (isoprene units)
    isoprene_pattern = Chem.MolFromSmarts("[C@H]([C@H](C)CC=C)")
    matches = mol.GetSubstructMatches(isoprene_pattern)
    if not matches:
        return False, "No terpenoid-like isoprene units found"

    # Check for rearranged/modified tetraterpene skeleton
    largest_aliphatic_ring = max(len(ring) for ring in mol.GetRingInfo().AtomRings())
    if largest_aliphatic_ring > 6:
        return True, "Contains a rearranged or modified tetraterpene skeleton"

    # If all checks pass, classify as a tetraterpenoid
    return True, "Contains a C40 skeleton with terpenoid-like isoprene units"