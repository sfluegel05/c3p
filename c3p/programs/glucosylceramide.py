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

    # Identify β-D-glucose moiety
    glucose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@H](CO)O1")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No β-D-glucose moiety found"

    # Identify amide linkage (N-C(=O))
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    # Identify sphingosine backbone: C18 chain with amine, two hydroxyls, and possible double bond
    sphingosine_patterns = [
        Chem.MolFromSmarts("[NH2]CC(CO)[C@H](O)CCCCCCCC=CCCC"),  # Common sphingosine pattern
        Chem.MolFromSmarts("[NH2]CC(CO)[C@H](O)CCCCCCCCCCCCCC") # Stereochemistry variant
    ]
    if not any(mol.HasSubstructMatch(pat) for pat in sphingosine_patterns):
        return False, "No sphingosine backbone found"

    # Identify long fatty acid chain (At least 16 carbons, preferably attached to the amide)
    long_chain_pattern = Chem.MolFromSmarts("C" * 16)  # Flexible long chain
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No sufficiently long fatty acid chain found"

    return True, "Contains glucose moiety linked to sphingosine backbone with a fatty acid chain"