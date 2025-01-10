"""
Classifies: CHEBI:52221 isothiocyanate
"""
from rdkit import Chem


def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    Isothiocyanates have the functional group N=C=S.

    Additional context around this group is checked to avoid false positives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isothiocyanate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the isothiocyanate pattern (N=C=S)
    isothiocyanate_pattern = Chem.MolFromSmarts("[$(C=N=C=S),$(N=C=S),$(S=C=N)]")

    # Check for the isothiocyanate functional group in a valid configuration
    if mol.HasSubstructMatch(isothiocyanate_pattern):
        # Ensure the group is properly bonded in a straightforward context
        atom_matches = mol.GetSubstructMatches(isothiocyanate_pattern)
        for match in atom_matches:
            atoms = [mol.GetAtomWithIdx(idx) for idx in match]
            nitrogen, carbon, sulfur  = atoms

            # Check if carbon has more than two other attachments (other than N=C=S itself)
            if carbon.GetDegree() > 2:
                return False, "N=C=S exists but context is chemically non-standard"
        
        return True, "Contains isothiocyanate functional group (N=C=S) in proper context"
    else:
        return False, "Does not contain isothiocyanate functional group (N=C=S)"