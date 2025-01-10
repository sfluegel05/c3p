"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate consists of a polyprenol chain linked to a phosphate or diphosphate group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise.
        str: Reason for classification or failure.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded isoprene unit pattern to recognize both cis and trans configurations
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C=CCC")
    
    # Check for repeated isoprene units
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 4:  # Improved heuristic for polyprenol length
        return False, "Insufficient isoprene units for a polyprenol chain"

    # Pattern for phosphate group (including variations of diphosphate linkage)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    diphosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-])OP(=O)([O-])[O-]")
    
    # Ensure phosphate (or diphosphate) group present
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    if not phosphate_matches and not diphosphate_matches:
        return False, "No phosphate or diphosphate group found"

    # Ensure phosphate or diphosphate is structurally positioned at the molecule end
    end_atoms = {match[0] for match in isoprene_matches}.union({match[-1] for match in isoprene_matches})
    attached_phosphate = False

    for match in phosphate_matches + diphosphate_matches:
        if set(match).intersection(end_atoms):
            attached_phosphate = True
            break

    if attached_phosphate:
        return True, "Valid polyprenol phosphate classified"

    return False, "Phosphate group does not attach correctly to the polyprenol chain"