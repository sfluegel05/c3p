"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: CHEBI:27849 iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdFMCS

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
    
    # Check for monoterpenoid characteristics
    n_carbon = rdMolDescriptors.CalcNumAtoms(mol, onlyObayedHeteroatoms=True, onlyHeavy=True, heteroatoms=[6])
    n_oxygen = rdMolDescriptors.CalcNumAtoms(mol, onlyObayedHeteroatoms=True, onlyHeavy=True, heteroatoms=[8])
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    
    if n_carbon < 10 or n_carbon > 15:
        return False, "Number of carbons outside typical range for monoterpenoids"
    if n_oxygen < 1 or n_oxygen > 3:
        return False, "Number of oxygens outside typical range for iridoid monoterpenoids"
    if n_rings < 2 or n_rings > 4:
        return False, "Number of rings outside typical range for iridoid monoterpenoids"
    
    # Check for additional structural features
    has_cyclopentane = mol.HasSubstructMatch(Chem.MolFromSmarts("[C@]1([C@@H](C)C)[C@@H](C)[C@@H](C)[C@@H]1"))
    has_pyran = mol.HasSubstructMatch(Chem.MolFromSmarts("C1=CCOCC1"))
    has_ester = mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)OC"))
    
    if not has_cyclopentane:
        return False, "No cyclopentane ring found"
    if not has_pyran:
        return False, "No pyran ring found"
    
    # Calculate confidence score
    confidence = 0
    if has_ester:
        confidence += 0.2
    if n_oxygen == 2:
        confidence += 0.3
    if n_rings == 3:
        confidence += 0.3
    
    mcs = rdFMCS.FindMCS([mol, iridoid_core], matchValences=False, ringMatchesRingOnly=True)
    mcs_ratio = mcs.numBonds / (mol.GetNumBonds() + iridoid_core.GetNumBonds() - mcs.numBonds)
    confidence += mcs_ratio * 0.2
    
    if confidence > 0.8:
        return True, f"Iridoid monoterpenoid with confidence score {confidence:.2f}"
    else:
        return False, f"Not a confident iridoid monoterpenoid (confidence {confidence:.2f})"