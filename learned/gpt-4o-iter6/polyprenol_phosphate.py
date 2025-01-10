"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate consists of a polyprenol chain linked to a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification or failure
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for isoprene unit (C=C-C-C=C)
    isoprene_pattern = Chem.MolFromSmarts("[CX3]=[CX3][CX3][CX3]=[CX3]")
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No polyprenol chain found, missing isoprene unit pattern"

    # Pattern for phosphoric acid ester, e.g., O-P(=O)(O)-
    phosphate_pattern = Chem.MolFromSmarts("O-P(=O)(O)-")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate ester group found"

    # Ensure phosphate is at the end of polyprenol (though positional specific check might be complex)
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)

    # Simple positional check: ensure a phosphate group is on any terminal end
    # Note: This assumes linear topology and might not handle cyclic variations
    ends = [0, len(mol.GetAtoms()) - 1]
    for end in ends:
        if any(end in match for match in phosphate_matches):
            return True, "Contains polyprenol chain with attached phosphate group"

    return False, "Phosphate group does not attach correctly to the polyprenol chain"