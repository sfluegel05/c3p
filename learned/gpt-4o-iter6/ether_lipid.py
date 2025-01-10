"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid has a glycerol-like backbone with one or more carbon atoms bonded to alkyl chains via ether linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adjusted glycerol-like backbone pattern with potential ether linkages
    glycerol_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][OX2]")  # Simplified for ethers
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol-like backbone found"

    # General ether linkage pattern
    ether_pattern = Chem.MolFromSmarts("[CX4]O[CX4]")  # General to detect ether connections
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if not ether_matches:
        return False, "No ether linkages found"

    # Check for phosphate groups (optional)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)(O)(O)[O-]")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)

    if has_phosphate:
        return True, "Contains ether linkage with glycerol-like backbone and phosphate group"

    return True, "Contains ether linkage with glycerol-like backbone, but no phosphate group found"

# Example test
# result, reason = is_ether_lipid("CCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(C)=O")
# print(result, reason)