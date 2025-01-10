"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid is any terpenoid derived from a diterpene (C20 skeleton),
    which may be rearranged or modified by removal of skeletal atoms.

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

    # Check for presence of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 22:
        return False, f"Number of carbons ({c_count}) not typical for diterpenoids (18-22)"

    # Check for isoprene units (C5 units)
    # Isoprene unit SMARTS pattern: C=C-C-C=C
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C[C@]C(=C)")
    matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(matches) < 2:
        # Try a looser pattern matching for isoprene units
        isoprene_pattern_loose = Chem.MolFromSmarts("C=C-C-C")
        matches_loose = mol.GetSubstructMatches(isoprene_pattern_loose)
        if len(matches_loose) < 4:
            return False, "Insufficient isoprene units found"

    # Check for ring systems (common in diterpenoids)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 1:
        return False, "No ring structures found, diterpenoids usually contain rings"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 450:
        return False, f"Molecular weight ({mol_wt:.2f} Da) not in typical diterpenoid range (250-450 Da)"

    # Check for oxygen atoms (terpenoids are terpenes with oxygen-containing functional groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found, terpenoids typically contain oxygen"

    # If all checks pass, classify as diterpenoid
    return True, "Molecule matches characteristics of a diterpenoid"