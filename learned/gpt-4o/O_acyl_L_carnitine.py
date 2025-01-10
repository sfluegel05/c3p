"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # O-acyl-L-carnitine pattern: L-carnitine component with specific chiral configuration and ester linkage
    o_acyl_l_carnitine_pattern = Chem.MolFromSmarts("O[C@@H](CC(=O)[O-])C[N+](C)(C)C")
    if not mol.HasSubstructMatch(o_acyl_l_carnitine_pattern):
        return False, "Does not match O-acyl-L-carnitine pattern"

    # Verify presence of ester linkage
    ester_pattern = Chem.MolFromSmarts("O=C-O[C@@H]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Ensure chiral configuration matches L-carnitine
    chiral_match = mol.HasSubstructMatch(o_acyl_l_carnitine_pattern)
    if not chiral_match:
        return False, "Chiral configuration does not match L-carnitine"

    return True, "Contains characteristic O-acyl-L-carnitine pattern with L-chirality"