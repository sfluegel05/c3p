"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid has one or more carbon atoms in glycerol bonded to an alkyl chain via ether linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated glycerol backbone pattern to include stereochemistry and complex variations
    glycerol_pattern = Chem.MolFromSmarts("[O][C@H][CH2][O]")  # Example pattern for stereo and oxygen groups
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for ether linkage (C-O-[!C])
    ether_pattern = Chem.MolFromSmarts("[CX4][OX2H0][!#1]")  # ether linkage with non-hydrogen atom
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if not ether_matches:
        return False, "No ether linkages found"

    # Count ether vs ester groups (O=C-O)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Ensuring more or equivalent ethers to esters for classification
    if not (len(ether_matches) > len(ester_matches)):
        return False, "Not enough ether linkages compared to ester linkages"

    # Look for phosphate groups (P(=O)(O)O)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)(O)O")
    if mol.HasSubstructMatch(phosphate_pattern):
        return True, "Contains ether linkage with glycerol backbone and phosphate group"

    return True, "Contains ether linkage with glycerol backbone, but no phosphate group found"

# You can test the function with examples:
# result, reason = is_ether_lipid("CCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(C)=O")
# print(result, reason)  # Expected output: True, "Contains ether linkage with glycerol backbone and phosphate group"