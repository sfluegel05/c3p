"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    Ether lipids have at least one carbon atom on glycerol bonded to an alkyl chain via an ether linkage,
    and often possess specific headgroups like phosphate or choline.

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

    # Pattern for alkyl ether linkage specifically in ether lipids
    ether_pattern = Chem.MolFromSmarts("C-O-[CH2]-[CH]")

    # Check for ether linkage pattern
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) < 1:
        return False, "No appropriate ether linkage found"

    # Check for common lipid headgroups (e.g., phosphate or choline)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)[O-]")
    choline_pattern = Chem.MolFromSmarts("[N+](C)(C)C")
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    has_choline = mol.HasSubstructMatch(choline_pattern)
    has_ethanolamine = mol.HasSubstructMatch(ethanolamine_pattern)

    # Ensure at least one typical headgroup is attached
    if has_phosphate or has_choline or has_ethanolamine:
        return True, "Contains ether linkages with common lipid headgroups such as phosphate, choline, or ethanolamine"
    else:
        return False, "Missing common lipid headgroups"