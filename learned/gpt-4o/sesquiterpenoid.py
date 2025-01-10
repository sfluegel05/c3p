"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is derived from a sesquiterpene with structural modifications
    possible. Typically a C15 skeleton, but may involve rearrangements or minor
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
    
    # Check for carbon backbone size (flexibility in carbon count due to rearrangements)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (13 <= carbon_count <= 15):
        return False, f"Expected approximately C15 count, found {carbon_count}"

    # Count for strategic sesquiterpenoid features
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count < 1:
        return False, "Sesquiterpenoids typically feature rings"

    # Check for typical functional groups (broadening scope)
    potential_groups = ['[OX2H]','[OX1C]([CX3]=[OX1])','[O]','C=O']  # Hydroxyls, carbonyls, epoxy
    group_found = any(mol.HasSubstructMatch(Chem.MolFromSmarts(group)) for group in potential_groups)
    if not group_found:
        return False, "No typical sesquiterpenoid functional groups found"
    
    # Check for complex arrangements (including simple rearrangements)
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if num_rotatable_bonds > 2:
        return True, "Possible sesquiterpenoid with characteristic chemical complexity"

    return True, "C15 backbone with sesquiterpenoid characteristics including possible functional or structural modifications"