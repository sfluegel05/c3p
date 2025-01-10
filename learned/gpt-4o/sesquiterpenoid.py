"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is derived from a sesquiterpene with structural modifications
    possible. Typically a C15 skeleton but may involve rearrangements or minor
    subtractive modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carbon backbone size with flexible range
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (12 <= carbon_count <= 16):  # Allowing for rearrangements
        return False, f"Expected approximately C15 count, found {carbon_count}"

    # Check for minimum number of rings to suggest cyclization
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count < 1:
        return False, "Typical sesquiterpenoids feature rings"
    
    # Detect key sesquiterpenoid functional groups using SMARTS patterns
    functional_groups = ['[OX2H]', '[C=O]', '[OX1C](=[OX1])', '[C-0]!@[C-0]']
    if not any(mol.HasSubstructMatch(Chem.MolFromSmarts(group)) for group in functional_groups):
        return False, "No typical sesquiterpenoid functional groups found"

    # Look for characteristic chemical complexity
    num_stereocenters = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    if num_stereocenters < 3:
        return False, "Low number of stereocenters, uncharacteristic of sesquiterpenoids"
    
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if num_rotatable_bonds > 3:
        return True, "Possible sesquiterpenoid with characteristic chemical complexity"

    return True, "C15-like backbone with sesquiterpenoid characteristics including functional or structural modifications"