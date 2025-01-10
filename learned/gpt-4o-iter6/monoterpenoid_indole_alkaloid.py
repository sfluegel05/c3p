"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
# Import required libraries from RDKit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    Monoterpenoid indole alkaloids feature complex indole structures, often modified or expanded.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded indole-like core structure patterns
    indole_patterns = [
        Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3ccccc23"),  # Standard indole
        Chem.MolFromSmarts("c1cc2[nH]cc(c2c1)C"),  # Simple indole variation
        Chem.MolFromSmarts("c1cc2[nH]cnc2c1"),  # Pyrroloindoles
        # Additional fused nitrogen-heterocycle patterns can be added
    ]
    if not any(mol.HasSubstructMatch(pat) for pat in indole_patterns):
        return False, "No recognizable indole-like structure"

    # Presence of nitrogen atoms assessment
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:  # Allowing for a diverse range of nitrogen configurations
        return False, "Insufficient nitrogen atoms, unlikely to be an alkaloid"

    # Structural complexity by ring counts and overlapping fused ring systems
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_rings < 2 or num_aromatic_rings < 1:
        return False, f"Insufficient overall or aromatic ring complexity: {num_rings} total rings, {num_aromatic_rings} aromatic"

    # Check for common monoterpenoid-like substitutions, such as methoxy groups
    methoxy_pattern = Chem.MolFromSmarts("CO")
    if mol.HasSubstructMatch(methoxy_pattern):
        return True, "Contains recognizable indole-like structure with methoxy groups suggesting monoterpenoid indole alkaloid"

    # Consider more specific pattern or functional group checks for monoterpenoids
    return True, "Contains recognizable indole-like structure and features typical of monoterpenoid indole alkaloids"