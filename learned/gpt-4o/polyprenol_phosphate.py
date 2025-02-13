"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate has a polyprenol chain attached via a phosphate ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group pattern
    phosphate_pattern = Chem.MolFromSmarts("P(O)(O)=O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for ester linkage in phosphate
    ester_linkage_pattern = Chem.MolFromSmarts("O-P")
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage to phosphate group found"

    # Check for the isoprene units (e.g., -C=C-CCC- pattern repeated)
    isoprene_pattern = Chem.MolFromSmarts("C=C-CC")
    isoprene_count = len(mol.GetSubstructMatches(isoprene_pattern))
    if isoprene_count < 2:  # Typically polyprenols have multiple repeating isoprene units
        return False, f"Insufficient isoprene units: {isoprene_count}"

    return True, "Contains polyprenol chain with phosphate ester bond"

# Note: This might miss some more complex polyprenol phosphates or give false results for similar structures