"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester has the general structure R-O-C(=O)-O-R', where R and R' are organyl groups.
    This also covers cyclic carbonates.

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

    # Define the core carbonate group SMARTS pattern where oxygens are connected to other atoms (R-O-C(=O)-O-R')
    carbonate_pattern = Chem.MolFromSmarts("[O](C=O)[O]")

    # Check if the molecule has the carbonate group
    if not mol.HasSubstructMatch(carbonate_pattern):
        return False, "No carbonate group found"

    # Get the substructure matches
    matches = mol.GetSubstructMatches(carbonate_pattern)
    for match in matches:
        carbon_idx = match[1]  # index of central carbon
        for oxygen_idx in [match[0], match[2]]: # indices of the two oxygens
            oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
            # check the atoms bonded to the oxygens of the carbonate. These must NOT be hydrogen atoms
            for bond in oxygen_atom.GetBonds():
                linked_atom = bond.GetOtherAtom(oxygen_atom)
                if linked_atom.GetAtomicNum() == 1:
                   return False, "One or more O groups linked to a Hydrogen"
    
    return True, "Contains a carbonate group with organyl groups attached."