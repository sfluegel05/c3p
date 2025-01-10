"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is a prenol phosphate resulting from the formal condensation 
    of the terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.

    Args:
        smiles (str): SMILES string

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define phosphate group pattern (monophosphate only)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)[O-0]")
    if phosphate_pattern is None:
        return False, "Invalid phosphate SMARTS pattern"

    # Search for phosphate group
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Define isoprene unit pattern: C=C-C-C (with appropriate atom types)
    isoprene_smarts = "[CH2]=[CH][CH2][CH2]"
    isoprene_pattern = Chem.MolFromSmarts(isoprene_smarts)
    if isoprene_pattern is None:
        return False, "Invalid isoprene SMARTS pattern"

    # Count isoprene units
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    num_isoprene_units = len(isoprene_matches)

    if num_isoprene_units >= 4:
        return True, f"Molecule is a polyprenol phosphate with {num_isoprene_units} isoprene units"
    else:
        return False, f"Found {num_isoprene_units} isoprene units, which is less than required"