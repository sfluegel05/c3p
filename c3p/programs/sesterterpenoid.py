"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid is derived from a sesterterpene, often featuring rearranged or modified C25 backbones.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Loosen the carbon count filter to range from 20 to 50
    if c_count < 20 or c_count > 50:
        return False, f"Carbon count of {c_count} is outside typical range for sesterterpenoids"

    # Check for possible terpenoid structural patterns; more generic recognition
    terpenoid_patterns = [
        Chem.MolFromSmarts("[CX3]=[CX3]"),  # Alternative double bond pattern often seen in terpenoids
        Chem.MolFromSmarts("C1=CC=C1"),     # Aromatic ring, as rearrangements might form such structures
    ]
    
    has_terpenoid_signature = any(mol.HasSubstructMatch(p) for p in terpenoid_patterns)
    
    if not has_terpenoid_signature:
        return False, "Does not have clear terpenoid-like structural features"

    # Acknowledge diversity in functional groups; hence, not strictly classifying based on specific groups
    return True, "Contains characteristics typical of a sesterterpenoid with a modified and flexible outlook on carbon skeleton and functional diversity"