"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyphenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for multiple aromatic rings containing hydroxyl groups
    phenol_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")  # Phenol or aromatic hydroxyl group
    matches = mol.GetSubstructMatches(phenol_pattern)
    if len(matches) < 2:
        return False, f"Found {len(matches)} phenolic rings, need at least 2 for polyphenol"

    # Further confirm by counting the total number of aromatic rings
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 2:
        return False, f"Found {num_aromatic_rings} aromatic rings, need at least 2 for polyphenol"

    return True, "Contains multiple phenolic rings, consistent with a polyphenol"