"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: CHEBI:33859 carbamate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester contains the pattern -O-C(=O)-N- where the nitrogen can be substituted.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carbamate pattern: -O-C(=O)-N-
    # [OX2] : oxygen with 2 connections
    # [CX3] : carbon with 3 connections
    # [OX1] : oxygen with 1 connection (double bonded)
    # [NX3] : nitrogen with 3 connections
    carbamate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[NX3]")
    
    # Alternative pattern for primary carbamates (-NH2)
    primary_carbamate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[NX2H2]")
    
    matches = mol.GetSubstructMatches(carbamate_pattern)
    primary_matches = mol.GetSubstructMatches(primary_carbamate_pattern)
    
    total_matches = len(matches) + len(primary_matches)
    
    if total_matches == 0:
        return False, "No carbamate ester group found"

    # Make sure the oxygen is not part of a carbonate (-O-C(=O)-O-) 
    carbonate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[OX2]")
    carbonate_matches = mol.GetSubstructMatches(carbonate_pattern)
    
    # For each carbamate match, verify it's not part of a carbonate
    valid_carbamates = 0
    for match in matches + primary_matches:
        is_carbonate = False
        for carbonate in carbonate_matches:
            if match[0] == carbonate[0] and match[1] == carbonate[1]:
                is_carbonate = True
                break
        if not is_carbonate:
            valid_carbamates += 1

    if valid_carbamates == 0:
        return False, "Found matching pattern but it's part of a carbonate group"

    # Check that the nitrogen is not part of an N-oxide
    n_oxide_pattern = Chem.MolFromSmarts("[NX4+]")
    if mol.HasSubstructMatch(n_oxide_pattern):
        for match in matches + primary_matches:
            n_atom = mol.GetAtomWithIdx(match[3])
            if n_atom.GetFormalCharge() > 0:
                valid_carbamates -= 1

    if valid_carbamates == 0:
        return False, "Found matching pattern but nitrogen is part of N-oxide"

    return True, f"Contains {valid_carbamates} carbamate ester group(s)"