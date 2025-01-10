"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    Tocols have a chromanol core with a hydrocarbon chain attached at position 2
    consisting of three isoprenoid units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved chroman-6-ol core pattern with stereochemistry consideration
    chromanol_pattern = Chem.MolFromSmarts("O1CC[C@](C2=CC=C(C=C2)C)([C@H]1)O")
    if not mol.HasSubstructMatch(chromanol_pattern):
        return False, "No chroman-6-ol core found"

    # Improved isoprenoid unit pattern; considering variability in saturation
    isoprenoid_unit = Chem.MolFromSmarts("C(=C)(C)CC")
    matches = mol.GetSubstructMatches(isoprenoid_unit)
    if len(matches) < 3:
        return False, f"Found {len(matches)} isoprenoid units, need at least 3"

    # Explicit verification of correct attachment at position 2
    chroman_6_ol_matches = mol.GetSubstructMatch(chromanol_pattern)
    if not chroman_6_ol_matches:
        return False, "Chromanol core detected but verification failed"

    # Ensure the chain spans three isoprenoid units
    carbon_2_index = chroman_6_ol_matches[2]  # Assuming index 2 corresponds to position 2
    carbon_2_atom = mol.GetAtomWithIdx(carbon_2_index)
    neighbors = [neighbor for neighbor in carbon_2_atom.GetNeighbors() if neighbor.GetSymbol() == 'C']
    
    # Check if one of the neighboring carbon atoms connects to the expected chain length
    valid_chain = False
    for neighbor in neighbors:
        path = AllChem.GetShortestPath(mol, carbon_2_index, neighbor.GetIdx())
        chain_length = len(path) - 1
        if chain_length >= 15:  # Minimum chain length for three isoprenoid units
            valid_chain = True
            break

    if not valid_chain:
        return False, "No valid hydrocarbon chain substitution at position 2 spanning three isoprenoid units"

    return True, "Contains chroman-6-ol core with valid hydrocarbon chain substitution spanning three isoprenoid units"