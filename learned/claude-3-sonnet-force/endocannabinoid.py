"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: CHEBI:37441 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids are cannabinoids present in mammalian biological fluids and tissues
    that activate cannabinoid receptors. They typically consist of a long fatty acid chain
    with an ethanolamine or glycerol head group, and may contain additional functional groups
    like epoxide rings or aromatic systems.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for fatty acid chain (long carbon chain with carboxyl group)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) != 1:
        return False, "No fatty acid chain found"

    # Look for ethanolamine or glycerol head group
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    glycerol_pattern = Chem.MolFromSmarts("[OX2][CX4]([CHX4])[OX2][CHX4][OX2]")
    head_group_matches = mol.GetSubstructMatches(ethanolamine_pattern) + mol.GetSubstructMatches(glycerol_pattern)
    if len(head_group_matches) != 1:
        return False, "No ethanolamine or glycerol head group found"

    # Look for epoxide ring
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)

    # Look for aromatic system
    aromatic_pattern = Chem.MolFromSmarts("a")
    aromatic_matches = mol.GetSubstructMatches(aromatic_pattern)

    # Check for long carbon chain (at least 16 carbons)
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 16:
        return False, "Carbon chain too short for endocannabinoid"

    # Classify based on structural features
    if len(epoxide_matches) > 0:
        return True, "Endocannabinoid with epoxide ring"
    elif len(aromatic_matches) > 0:
        return True, "Endocannabinoid with aromatic system"
    else:
        return True, "Endocannabinoid with long fatty acid chain"