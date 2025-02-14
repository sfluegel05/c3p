"""
Classifies: CHEBI:17297 UDP-sugar
"""
"""
Classifies: CHEBI:27599 UDP-sugar

A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to an unspecified sugar via an anomeric diphosphate linkage.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for pyrimidine ring pattern
    pyrimidine_pattern = Chem.MolFromSmarts("nc1ncnc(n1)")
    pyrimidine_match = mol.GetSubstructMatches(pyrimidine_pattern)
    if not pyrimidine_match:
        return False, "No pyrimidine ring found"
    
    # Look for ribose/deoxyribose sugar pattern
    sugar_pattern = Chem.MolFromSmarts("[OX2r]~[CX4r]~[OX2r]~[CX4r]~[OX2r]")
    sugar_match = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_match:
        return False, "No sugar substructure found"
    
    # Look for diphosphate linkage pattern
    diphosphate_pattern = Chem.MolFromSmarts("[OX2r][P]([OX2r])([OX2r])=O")
    diphosphate_match = mol.GetSubstructMatches(diphosphate_pattern)
    if not diphosphate_match:
        return False, "No diphosphate linkage found"
    
    # Check if all components (pyrimidine, sugar, diphosphate) are present
    if pyrimidine_match and sugar_match and diphosphate_match:
        return True, "Contains a pyrimidine, a sugar, and a diphosphate linkage (UDP-sugar)"
    
    return False, "Missing one or more components of a UDP-sugar"