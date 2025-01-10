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

    # Define improved isoprene unit patterns. These include more flexible configurations.
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C=C")  # Simplified for linking head-to-tail
    # Check for repeated isoprene units
    if mol.GetNumAtoms() < 15:  # Assuming polyprenols are larger molecules with multiple isoprene units
        return False, "Molecule too small to contain a significant polyprenol chain"

    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No significant polyprenol chain found, missing isoprene repeats"

    # Pattern for phosphoric acid ester, e.g., phosphate group
    phosphate_pattern = Chem.MolFromSmarts("O-P(=O)(O)-")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate ester group found"

    # Ensure phosphate is structurally positioned at the molecule end
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    num_isoprene_units = len(mol.GetSubstructMatches(isoprene_pattern))

    # Neglect molecules with fewer than a polynomial apt isoprene chain (heuristic)
    if num_isoprene_units < 3:
        return False, "Insufficient isoprene units for classification as a polyprenol"

    # Check for terminal phosphate connection (based localization and match, assuming it end-caps an isoprene chain setup)
    atom_count = len(mol.GetAtoms())  # crude size estimate
    end_atoms = [0, atom_count - 1]

    for end in end_atoms:
        if any(end in match for match in phosphate_matches):
            return True, "Valid polyprenol phosphate classified"

    return False, "Phosphate group does not attach correctly to the polyprenol chain"