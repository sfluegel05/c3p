"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA is characterized by a 3-hydroxy fatty acid linked to coenzyme A.

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

    # SMARTS pattern for Coenzyme A moiety
    coA_pattern = Chem.MolFromSmarts("[SH]C(=O)NCCNC(=O)C1=NC(N)=NC2=C1N=CN=C2N3C[C@H]4O[C@@H](COP(O)(=O)OP(O)(O)=O)[C@H]4O[C@H]3CO")

    # SMARTS pattern for 3-hydroxy fatty acid moiety
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("CC(O)C(=O)")

    # Check for Coenzyme A moiety
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "No Coenzyme A moiety found"

    # Check for 3-hydroxy fatty acid component
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "No 3-hydroxy fatty acid moiety found"

    return True, "Molecule contains both Coenzyme A and 3-hydroxy fatty acid moieties"