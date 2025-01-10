"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    Gangliosides are glycosphingolipids that contain one or more sialic acid residues
    linked to the sugar chain, and a ceramide backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined pattern for ceramide backbone with sphingosine structure
    ceramide_pattern = Chem.MolFromSmarts("C(=O)NC(CO)C\C=C\CCCCC")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No valid ceramide backbone found"

    # Enhanced pattern for sialic acid residues substructure
    sialic_acid_pattern = Chem.MolFromSmarts("C[C@H](O)C(=O)OC")  # Accommodates variations in sialic acids
    if not mol.HasSubstructMatch(sialic_acid_pattern):
        return False, "No valid sialic acid residue found"

    # Sufficient variety and linkage of sugar residues (flexible glycans)
    sugar_pattern_1 = Chem.MolFromSmarts("COC[C@H]1O[C@H](COP(O)(=O)O)[C@H](O)[C@@H]1O")
    sugar_pattern_2 = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1")  # Generic pyranose form
    sugar_matches_1 = mol.GetSubstructMatches(sugar_pattern_1)
    sugar_matches_2 = mol.GetSubstructMatches(sugar_pattern_2)
    if len(sugar_matches_1) + len(sugar_matches_2) < 3:
        return False, "Missing necessary glycan residues"

    return True, "Contains ceramide backbone with sialic acid and necessary glycan residues"