"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    Ether lipids have at least one carbon atom on glycerol bonded to an alkyl chain via an ether linkage,
    and often possess specific headgroups like phosphate, choline, or ethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized pattern for ether linkage in glycerolipids
    ether_linkage_patterns = [
        Chem.MolFromSmarts("C-O-[CH2]-[CH2]"),  # Simple ether linkage
        Chem.MolFromSmarts("[C;H2]-[O]-[C;H2]")  # More general ether pattern accounting for various configurations
    ]

    # Check for ether linkage pattern
    if not any(mol.HasSubstructMatch(pattern) for pattern in ether_linkage_patterns):
        return False, "No appropriate ether linkage found"

    # Check for common lipid headgroups patterns more broadly
    headgroup_patterns = [
        Chem.MolFromSmarts("P(=O)(O)[O-]"),  # Phosphate
        Chem.MolFromSmarts("[N+](C)(C)C"),   # Choline
        Chem.MolFromSmarts("NCCO"),          # Ethanolamine
        Chem.MolFromSmarts("NC(=O)C")        # Possible amide linkage
        # Add more headgroup patterns as necessary for the lipid class
    ]
    
    # Ensure the presence of at least one typical lipid headgroup
    if any(mol.HasSubstructMatch(pattern) for pattern in headgroup_patterns):
        return True, "Contains ether linkages with common lipid headgroups such as phosphate, choline, or ethanolamine"
    else:
        return False, "Missing common lipid headgroups"