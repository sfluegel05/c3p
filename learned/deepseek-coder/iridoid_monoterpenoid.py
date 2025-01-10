"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: CHEBI:23835 iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    Iridoid monoterpenoids are characterized by a cyclopentane ring fused to a six-membered oxygen heterocycle.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more specific SMARTS pattern for the fused cyclopentane and six-membered oxygen heterocycle
    iridoid_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]1@[C]2[C][C][C][C][O]2")
    if not mol.HasSubstructMatch(iridoid_pattern):
        return False, "No cyclopentane ring fused to a six-membered oxygen heterocycle found"

    # Check for monoterpenoid backbone (typically 10 carbons, derived from isoprene units)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8 or c_count > 30:
        return False, f"Carbon count ({c_count}) not consistent with monoterpenoid backbone"

    # Check for oxygen heterocycle (at least one oxygen in a ring)
    o_in_ring = any(atom.GetAtomicNum() == 8 and atom.IsInRing() for atom in mol.GetAtoms())
    if not o_in_ring:
        return False, "No oxygen atom in a ring (required for oxygen heterocycle)"

    # Check molecular weight (iridoids typically have MW between 150-800 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 800:
        return False, f"Molecular weight ({mol_wt:.2f} Da) not consistent with iridoid monoterpenoids"

    # Check for isoprene-like substructures
    isoprene_pattern = Chem.MolFromSmarts("[C]=[C]-[C]=[C]")
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No isoprene-like substructure found"

    return True, "Contains cyclopentane ring fused to a six-membered oxygen heterocycle, consistent with iridoid monoterpenoid structure"