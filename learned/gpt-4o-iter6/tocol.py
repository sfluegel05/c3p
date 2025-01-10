"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    A tocol is defined as having a chromanol core with a hydrocarbon chain
    consisting of three isoprenoid units attached at position 2.

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
    
    # Define chroman-6-ol core pattern
    chromanol_pattern = Chem.MolFromSmarts("c1(ccc2c(c1)CC(O2)C)O")
    if not mol.HasSubstructMatch(chromanol_pattern):
        return False, "No chroman-6-ol core found"

    # Define isoprenoid pattern as a flexible five-carbon chain structure
    isoprenoid_pattern = Chem.MolFromSmarts("C(C)=C-C")
    matches = mol.GetSubstructMatches(isoprenoid_pattern)
    if len(matches) < 3:
        return False, f"Found {len(matches)} isoprenoid units, need at least 3"

    # Check if the hydrocarbon chain is correctly attached to the position 2 carbon
    chroman_6_ol_match = mol.GetSubstructMatch(chromanol_pattern)
    if chroman_6_ol_match:
        carbon_2_index = chroman_6_ol_match[1]  # position 2 in chroman-6-ol
        carbon_2_atom = mol.GetAtomWithIdx(carbon_2_index)
        # Check for carbon attachment
        if not any(neighbor.GetSymbol() == 'C' for neighbor in carbon_2_atom.GetNeighbors()):
            return False, "No valid hydrocarbon chain substitution at position 2"

    # Return true if all checks passed
    return True, "Contains chroman-6-ol core with valid hydrocarbon chain"