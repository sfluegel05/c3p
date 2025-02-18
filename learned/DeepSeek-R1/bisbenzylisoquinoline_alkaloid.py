"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: CHEBI:???? bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    These compounds consist of two benzylisoquinoline units linked by ether bridges.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Define isoquinoline core pattern (benzene fused to pyridine)
    isoquinoline_pattern = Chem.MolFromSmarts("c1ccc2cnccc2c1")
    # Check for at least two isoquinoline units
    isoq_matches = mol.GetSubstructMatches(isoquinoline_pattern)
    if len(isoq_matches) < 2:
        return False, f"Found {len(isoq_matches)} isoquinoline units, need at least 2"

    # Look for ether bridges (O connecting two aromatic carbons)
    ether_pattern = Chem.MolFromSmarts("[c,a]-[O]-[c,a]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) < 2:
        return False, f"Found {len(ether_matches)} ether bridges, need at least 2"

    # Check for tertiary amines (common in these alkaloids)
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3;H0;R]")
    amine_matches = mol.GetSubstructMatches(tertiary_amine_pattern)
    if len(amine_matches) < 2:
        return False, f"Found {len(amine_matches)} tertiary amines, need at least 2"

    # Check for methylenedioxy groups (O-C-O)
    methylenedioxy_pattern = Chem.MolFromSmarts("O-C-O")
    methylenedioxy_matches = mol.GetSubstructMatches(methylenedioxy_pattern)
    if methylenedioxy_matches:
        return True, "Contains two isoquinoline units, ether bridges, and methylenedioxy groups"

    # Check molecular weight (typically >500 for dimers)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight {mol_wt} too low for dimer"

    return True, "Contains two isoquinoline units linked by ether bridges"