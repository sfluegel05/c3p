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

    # Identify the chromanol core (chroman-6-ol)
    chromanol_pattern = Chem.MolFromSmarts("O1CC(C)C2=C(C1)C=C(C=C2)C")
    if not mol.HasSubstructMatch(chromanol_pattern):
        return False, "No chroman-6-ol core found"

    # Find the isoprenoid-like hydrocarbon chain (e.g., -CCC=C)
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("CCC=C")
    matches = mol.GetSubstructMatches(hydrocarbon_chain_pattern)
    if len(matches) < 3:
        return False, f"Found {len(matches)} isoprenoid units, need at least 3"

    # Check for the correct substitution at position 2
    # Assuming position 2 is the second carbon in the chromanol core:
    chain_attachment_site = chromanol_pattern.GetSubstructMatch(mol)[1]
    num_attached = sum(1 for atom in mol.GetAtomWithIdx(chain_attachment_site).GetNeighbors() if atom.GetSymbol() == 'C')
    if num_attached < 1:
        return False, "No valid hydrocarbon chain substitution at position 2"

    # Return true if all checks passed
    return True, "Contains chroman-6-ol core with valid hydrocarbon chain"