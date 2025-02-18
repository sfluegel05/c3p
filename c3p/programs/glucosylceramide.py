"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide contains a glucose moiety, a sphingosine backbone, and a long fatty acid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for β-D-glucose moiety using revised pattern
    glucose_pattern = Chem.MolFromSmarts("C(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No β-D-glucose moiety found"

    # Check for amide linkage attached to a long chain (not precise but helps identify linkage point)
    amide_pattern = Chem.MolFromSmarts("NC(=O)CCCCCCCCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage to a long fatty acid chain found"

    # Check for sphingosine backbone pattern (N-alkylated amine with long chain and specific OH groups)
    sphingosine_patterns = [
        Chem.MolFromSmarts("[NX3]C[C@H](O)CO"), # Basic pattern to identify the common N-linkage
    ]
    if not any(mol.HasSubstructMatch(pat) for pat in sphingosine_patterns):
        return False, "No sphingosine backbone found"

    return True, "Contains glucose moiety linked to a sphingosine-like backbone with a fatty acid chain"