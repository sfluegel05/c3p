"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is derived from a sesquiterpene with structural modifications
    possible. Typically C15 skeleton, with possible rearrangements or minor
    subtractive modifications.
    
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
    
    # Check if molecule has C15 backbone typically indicative of sesquiterpenoids
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 13 or carbon_count > 15:
        return False, f"Expected 13-15 carbons, found {carbon_count}"

    # Look for typical sesquiterpenoid features such as ring and possible terpenoid patterns
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count == 0:
        return False, "No rings present; unlikely to be a sesquiterpenoid"

    # Check for strategic double bonds or hydroxyl groups that are common
    hydroxy_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0)
    if hydroxy_count < 1:
        return False, "No hydroxyl groups found; sesquiterpenoids often have them"

    # Simple check for complex arrangements (could be extended by more detailed SMARTS patterns)
    possible_rearrangement = rdMolDescriptors.CalcNumRotatableBonds(mol) > 3
    
    if possible_rearrangement:
        return True, "C15 backbone with potential rearrangement features detected typical of sesquiterpenoids"

    return True, "C15 backbone with sesquiterpenoid characteristics"