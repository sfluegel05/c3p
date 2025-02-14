"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol phosphorylated at one or more hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a glycerol backbone related to phosphatidylinositol
    # Example pattern: C(CO)COC(=O)
    glycerol_pattern = Chem.MolFromSmarts("CC(O)COC(=O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone part of phosphatidylinositol found"

    # Check for inositol ring: A six-membered ring with multiple hydroxy groups
    # The stereo/ordering might vary
    inositol_ring_pattern = Chem.MolFromSmarts("C1(C(C(C(C(C1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(possible_inositol):
        return False, "No inositol ring with hydroxy groups found"

    # Check for the presence of phosphate groups attached to the inositol.
    # More than one phosphate group should be attached to different positions on the inositol.
    phosphate_group_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_group_pattern)

    # Check if phosphate groups are connected specifically to the inositol ring
    # This is tricky without labeled hydrogens, but checking adjacency may suffice
    adjacency_match = any(
        mol.GetAtomWithIdx(pm[0]).GetImplicitValence() > 1 for pm in phosphate_matches
    )
    if not adjacency_match:
        return False, "No phosphate groups properly attached to inositol"

    return True, "Contains glycerol backbone, inositol ring, and phosphorylated groups specific to phosphoinositides"