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
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized "glycerol-like" backbone pattern to include stereochemistry and additional groups
    glycerol_pattern = Chem.MolFromSmarts("[O][CH][CH2O][!#1]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol-like backbone found"

    # Revised ether linkage pattern suitable for complex structures
    ether_pattern = Chem.MolFromSmarts("[CH2]O[CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if not ether_matches:
        return False, "No ether linkages found"

    # Identifying ester groups which should be throughout molecule
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0][CX4]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Ensuring there are enough ethers to suggest classification
    if len(ether_matches) < 1:
        return False, "Not enough ether linkages present"

    # Determine if molecule contains a phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)(O)[O-]")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)

    if has_phosphate:
        return True, "Contains ether linkage with glycerol-like backbone and phosphate group"
    
    return True, "Contains ether linkage with glycerol-like backbone, but no phosphate group found"

# Test the function with examples
# result, reason = is_ether_lipid("CCCCCCCCCCCOC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(C)=O")
# print(result, reason)