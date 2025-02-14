"""
Classifies: CHEBI:37143 organofluorine compound
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound is defined as a compound containing at least one carbon-fluorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carbon-fluorine bond
    cf_pattern = Chem.MolFromSmarts("[C]F")
    
    if mol.HasSubstructMatch(cf_pattern):
        # Check if F is connected to a C
        matches = mol.GetSubstructMatches(cf_pattern)
        for match in matches:
            f_atom_idx = match[1]
            f_atom = mol.GetAtomWithIdx(f_atom_idx)
            for neighbor in f_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    return True, "Contains at least one carbon-fluorine bond"

    return False, "Does not contain a carbon-fluorine bond"