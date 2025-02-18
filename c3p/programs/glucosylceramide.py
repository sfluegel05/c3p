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

    # Check for β-D-glucose moiety using a more precise SMARTS pattern
    glucose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H](CO)[C@H](O)[C@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(glucose_pattern):
        return (False, "No β-D-glucose moiety found")

    # Check for an amide bond to a long fatty acid chain
    amide_pattern = Chem.MolFromSmarts("NC(=O)C[C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0]")
    if not mol.HasSubstructMatch(amide_pattern):
        return (False, "No amide linkage to a long fatty acid chain found")

    # Check for a basic sphingosine backbone pattern
    # We focus on characteristic features: C-OH groups with amine linkages, possibly with C=C for unsaturation
    sphingosine_pattern = Chem.MolFromSmarts("N[C;$(CO)][CX4H2][CX4H][CX4](OC)[CX4H2]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return (False, "No sphingosine backbone found")

    return (True, "Contains glucose moiety linked to a sphingosine-like backbone with a fatty acid chain")