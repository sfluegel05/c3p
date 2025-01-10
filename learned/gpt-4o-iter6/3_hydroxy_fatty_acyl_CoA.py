"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA consists of a 3-hydroxy fatty acid moiety linked to coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more comprehensive SMARTS pattern for CoA
    coA_pattern = Chem.MolFromSmarts("NC1=NC=NC2=C1N=C(NC2=O)N3[C@H]4O[C@H]([C@@H](O[C@H]4COP([O-])(=O)OP(O)(=O)O)O)C34O[C@@H]([C@H]4[OP(O)(O)=O])CO3")

    # Define a general pattern for a 3-hydroxy fatty acid
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("C[C@@H](O)CC(=O)")

    # Check for CoA moiety
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "No coenzyme A moiety found"

    # Check for 3-hydroxy fatty acid component
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "No 3-hydroxy fatty acid moiety found"

    return True, "Molecule contains both coenzyme A and 3-hydroxy fatty acid moieties"