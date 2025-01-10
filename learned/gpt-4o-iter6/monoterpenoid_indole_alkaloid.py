"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
# Import required libraries from RDKit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    Monoterpenoid indole alkaloids typically feature an indole ring fused with monoterpene-derived structures.

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

    # Define a SMARTS pattern to identify the indole core structure
    indole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3ccccc23")
    if not mol.HasSubstructMatch(indole_pattern):
        # Check for other fused or bridged indole variations as fallback
        alternative_indole_patterns = [
            Chem.MolFromSmarts("c1cc2[nH]cc(c2c1)C")  # Example variation: check for a simple indolic variation
            # Add more pattern variations as needed to cover known substructures
        ]
        if not any(mol.HasSubstructMatch(pat) for pat in alternative_indole_patterns):
            return False, "No recognizable indole-like structure"

    # Check for nitrogen presence, a common element in alkaloids
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:  # Typically more than one nitrogen in complex structures
        return False, "Insufficient nitrogen atoms, unlikely to be an alkaloid"

    # Verify structural complexity through additional rings complexity facing typical of monoterpenoid systems
    num_fused_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_fused_rings < 3:  # An arbitrary choice, assuming complex alkaloidal structure should have multiple rings
        return False, f"Found {num_fused_rings} rings, not complex enough"

    # For compound-specific motifs in monoterpenoid indole alkaloids, extra checks could be added here
    # E.g., check for specific functional groups or special ring junctions frequently seen in monoterpenoids

    return True, "Contains recognizable indole-like structure and features typical of monoterpenoid indole alkaloids"