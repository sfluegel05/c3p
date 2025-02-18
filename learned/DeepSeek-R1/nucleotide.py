"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide consists of a nucleoside (sugar + base) with at least one phosphate group
    attached to the 3' or 5' hydroxyl group of the sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for nucleoside: sugar (ribose-like) connected to a heterocyclic base
    # Sugar pattern: O-C1-C-O-C-O (approximate furanose structure)
    # Base: aromatic nitrogen-containing ring (e.g., purine/pyrimidine)
    nucleoside_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@@H](O)[C@H](O1)n2cnc3c2ncnc3N")
    if not mol.HasSubstructMatch(nucleoside_pattern):
        return False, "No nucleoside moiety detected"

    # Check for phosphate group attached to 3' or 5' oxygen
    # Phosphate ester pattern: O-P(=O)(O)-O connected to sugar
    phosphate_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@@H](OP(=O)(O)O)[C@H](O1)n")
    if not mol.HasSubstructMatch(phosphate_pattern):
        # Check alternative positions (5' or 3')
        alt_phosphate = Chem.MolFromSmarts("[C@H]1[C@H](OP(=O)(O)O)[C@@H](O)[C@H](O1)n")
        if not mol.HasSubstructMatch(alt_phosphate):
            return False, "No phosphate group attached to 3' or 5' position"

    return True, "Nucleoside with phosphate group at 3' or 5' position"