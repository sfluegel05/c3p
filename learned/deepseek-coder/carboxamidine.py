"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:37671 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine has the structure RC(=NR)NR2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxamidine pattern: RC(=NR)NR2
    # More specific pattern to avoid false positives
    carboxamidine_pattern = Chem.MolFromSmarts("[CX3](=[NX2])[NX3H2,NX3H1,NX3H0]")
    
    # Check for the basic pattern
    if not mol.HasSubstructMatch(carboxamidine_pattern):
        return False, "No carboxamidine structure (RC(=NR)NR2) found"

    # Additional checks to reduce false positives
    # 1. Ensure the nitrogen atoms are not part of a guanidine group
    guanidine_pattern = Chem.MolFromSmarts("[NX3H2,NX3H1,NX3H0][CX3](=[NX2])[NX3H2,NX3H1,NX3H0]")
    if mol.HasSubstructMatch(guanidine_pattern):
        return False, "Contains guanidine structure, not carboxamidine"

    # 2. Ensure the central carbon has exactly one double bond to nitrogen
    central_carbon_pattern = Chem.MolFromSmarts("[CX3](=[NX2])")
    matches = mol.GetSubstructMatches(central_carbon_pattern)
    if len(matches) != 1:
        return False, "Multiple potential carboxamidine centers found"

    # 3. Verify the connectivity of the atoms
    # Get the matched atoms
    match = mol.GetSubstructMatch(carboxamidine_pattern)
    if len(match) < 3:
        return False, "Invalid carboxamidine structure"

    # Check that the central carbon is bonded to exactly 3 atoms
    central_carbon = match[0]
    if len(mol.GetAtomWithIdx(central_carbon).GetBonds()) != 3:
        return False, "Central carbon has incorrect connectivity"

    return True, "Contains the carboxamidine structure (RC(=NR)NR2)"