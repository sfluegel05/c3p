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
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine ring found"
    
    # Look for ribose/deoxyribose sugar pattern
    sugar_pattern = Chem.MolFromSmarts("[OX2r]~[CX4r]~[OX2r]~[CX4r]~[OX2r]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar substructure found"
    
    # Look for diphosphate linkage pattern
    diphosphate_pattern = Chem.MolFromSmarts("[OX2r][P]([OX2r])([OX2r])=O")
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    if not diphosphate_matches:
        return False, "No diphosphate linkage found"
    
    # Check if the diphosphate connects pyrimidine and sugar
    for match in diphosphate_matches:
        pyrimidine_atom = mol.GetAtomWithIdx(match[2])
        sugar_atom = mol.GetAtomWithIdx(match[0])
        if pyrimidine_atom.IsInRing() and sugar_atom.IsInRing():
            return True, "Contains a pyrimidine linked to a sugar via an anomeric diphosphate linkage"
    
    return False, "Pyrimidine and sugar not connected via a diphosphate linkage"