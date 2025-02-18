"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide is characterized by a β-D-glucose moiety, a sphingosine-like backbone,
    and an amide-linked long fatty acid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Refined SMARTS pattern for β-D-glucose moiety
    glucose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H](O[C@H]1O)CO")
    if not mol.HasSubstructMatch(glucose_pattern):
        return (False, "No β-D-glucose moiety found")

    # Refined SMARTS pattern for an amide bond to a long fatty acid chain
    amide_pattern = Chem.MolFromSmarts("C(=O)N[C@@H](C=O)C[C;R0][C;R0]")  # Broadened pattern for the fatty chain
    if not mol.HasSubstructMatch(amide_pattern):
        return (False, "No amide linkage to a long fatty acid chain found")

    # Refined pattern for a sphingosine-like backbone; allows flexibility and common functional groups
    sphingosine_pattern = Chem.MolFromSmarts("[N;R0][CX4]-[CX4]-[CX4](-[CX4])[O]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return (False, "No sphingosine backbone found")

    return (True, "Contains glucose moiety linked to a sphingosine-like backbone with a fatty acid chain")