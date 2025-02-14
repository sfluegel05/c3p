"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester has the general structure R-O-C(=O)-O-R', where R and R' are organyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core carbonate group SMARTS pattern
    carbonate_pattern = Chem.MolFromSmarts("O=C(O)O")

    # Check if the molecule has the carbonate group.
    if not mol.HasSubstructMatch(carbonate_pattern):
        return False, "No carbonate group found"

    # Check that the groups attached to the oxygens are carbons, not hydrogens
    carbonate_matches = mol.GetSubstructMatches(carbonate_pattern)
    for match in carbonate_matches:
        carbon_idx = match[1] # The central carbon index
        for bond in mol.GetAtomWithIdx(carbon_idx).GetBonds():
            other_atom = bond.GetOtherAtom(mol.GetAtomWithIdx(carbon_idx))
            if other_atom.GetAtomicNum() == 8: # Check if linked to O
              for linked_bond in other_atom.GetBonds():
                linked_atom = linked_bond.GetOtherAtom(other_atom)
                if linked_atom.GetAtomicNum() != 1 and linked_atom.GetAtomicNum() != 6:  # Check that it is not hydrogen nor carbon
                  return False, "One or more O groups linked to non-carbon atoms"

    return True, "Contains a carbonate group with organyl groups attached."