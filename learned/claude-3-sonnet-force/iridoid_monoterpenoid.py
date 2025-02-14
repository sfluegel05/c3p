"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: CHEBI:27849 iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    Iridoid monoterpenoids are monoterpenoids biosynthesized from isoprene and often intermediates
    in the biosynthesis of alkaloids. They usually consist of a cyclopentane ring fused to a six-membered
    oxygen heterocycle; cleavage of a bond in the cyclopentane ring gives rise to the secoiridoids.

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

    # Normalize molecule
    mol = Chem.RemoveHs(mol)

    # Check for iridoid core structure
    iridoid_core = Chem.MolFromSmarts("[C@]12[C@@](C)(C[C@@H]1C)[C@@H](O)[C@@H](C2)O")
    if not mol.HasSubstructMatch(iridoid_core):
        return False, "No iridoid core structure found"

    # Check for glucose or glucopyranoside moiety
    glucose_pattern = Chem.MolFromSmarts("[C@H]1O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    has_glucose = mol.HasSubstructMatch(glucose_pattern)

    # Check for ester or acyl groups
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    acyl_pattern = Chem.MolFromSmarts("C(=O)C")
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_acyl = mol.HasSubstructMatch(acyl_pattern)

    # Check for molecular size and complexity
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol)

    if mol_weight < 300 or mol_weight > 700:
        return False, "Molecular weight outside typical range for iridoid monoterpenoids"
    if n_rings < 2 or n_rings > 5:
        return False, "Number of rings outside typical range for iridoid monoterpenoids"
    if n_rotatable < 3 or n_rotatable > 15:
        return False, "Number of rotatable bonds outside typical range for iridoid monoterpenoids"
    if n_heavy_atoms < 20 or n_heavy_atoms > 50:
        return False, "Number of heavy atoms outside typical range for iridoid monoterpenoids"

    # Calculate score
    score = 0
    if has_glucose:
        score += 0.3
    if has_ester or has_acyl:
        score += 0.2
    if n_rings == 3:
        score += 0.2
    if n_rotatable >= 8:
        score += 0.2
    mcs = rdFMCS.FindMCS([mol, iridoid_core], matchValences=False, ringMatchesRingOnly=True)
    mcs_ratio = mcs.numBonds / (mol.GetNumBonds() + iridoid_core.GetNumBonds() - mcs.numBonds)
    score += mcs_ratio * 0.1

    if score >= 0.8:
        return True, f"Iridoid monoterpenoid with score {score:.2f}"
    else:
        return False, f"Not a confident iridoid monoterpenoid (score {score:.2f})"