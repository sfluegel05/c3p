"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    A triterpenoid is characterized by a core C30 structure derived from a triterpene,
    which can include rearrangements or modifications by removal of skeletal atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (20 <= c_count <= 50):  # Expanded carbon range to catch variations and modifications
        return False, f"Carbon count {c_count} is not in the expected range for triterpenoids"

    # Allow flexibility in ring structures while expecting complex polycyclic arrangements
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2:  # Assuming at least two rings to cover some rearranged structures
        return False, "Too few rings for typical triterpenoid structures"

    # Include additional functional groups typical to triterpenoids
    functional_groups = [
        Chem.MolFromSmarts("C=O"),  # Carbonyl group
        Chem.MolFromSmarts("O"),    # Hydroxyl/methoxy/ether groups
        Chem.MolFromSmarts("CO"),   # Carboxylate ester structures
        Chem.MolFromSmarts("O=C(O)"), # Carboxylic acid group
        Chem.MolFromSmarts("C=C")   # Double bonds, indicative of unsaturation
    ]
    if not any(mol.HasSubstructMatch(fg) for fg in functional_groups):
        return False, "No characteristic triterpenoid functional groups found"

    return True, "Structure conforms with general criteria for triterpenoids"