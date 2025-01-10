"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is derived from a sesquiterpene, typically a C15 skeleton
    but may involve structural modifications such as rearrangements or subtractive
    modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for approximate carbon backbone size
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 10 or carbon_count > 25:  # Allow more flexibility
        return False, f"Expected more carbon atoms typical of sesquiterpenoids, found {carbon_count}"

    # Check for any ring structure
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count == 0:
        return False, "No ring structures, which are common in sesquiterpenoids"

    # Check for common sesquiterpenoid functionalities
    common_smarts = [
        '[C;X4]',  # Tetrahedral carbons (indicating complex alkyl chains)
        '[O][C]=O',  # Ester groups 
        '[O][#6]',   # Possible alcohols or ethers
        '[C](=O)[C]',  # Ketones often found in sesquiterpene lactones
    ]
    found_group = False
    for smarts in common_smarts:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            found_group = True
            break
    if not found_group:
        return False, "No typical sesquiterpenoid functional groups found"

    return True, "C15-like skeleton with sesquiterpenoid functional groups and possible structural modifications"