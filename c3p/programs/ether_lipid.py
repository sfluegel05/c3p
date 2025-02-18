"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    Ether lipids have one or more carbon atoms on glycerol bonded to an alkyl chain via an ether linkage.

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

    # Improved SMARTS pattern for glycerol backbone with flexible connectivity
    glycerol_pattern = Chem.MolFromSmarts("[CH2][O][CH][O][CH2]")

    # Adjust ether linkage patterns to cater for stereo and other forms of bonding
    ether_pattern = Chem.MolFromSmarts("[CX4,O]O[CX4]")

    # Look for ether linkage pattern
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) < 1:
        return False, "No ether linkage found"

    # Improved headgroup checks
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    choline_pattern = Chem.MolFromSmarts("C[N+](C)(C)C")
    
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    has_choline = mol.HasSubstructMatch(choline_pattern)

    if has_phosphate or has_choline:
        return True, "Contains ether linkages along with common lipid headgroups (e.g., phosphate, choline)"
    else:
        return False, "Missing common lipid headgroups"