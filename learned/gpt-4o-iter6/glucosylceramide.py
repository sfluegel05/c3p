"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide has a sphingosine backbone with an amide-linked fatty acyl chain
    and a β-D-glucose sugar attached to the primary hydroxyl group.

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

    # Define sphingosine backbone SMARTS with improved accuracy
    sphingosine_pattern = Chem.MolFromSmarts('N[C@@H](CO)C(O)C=CC')
    # Here, ensure patterns look for at least one double bond in the backbone
    # and check for the correct positions for hydroxy groups.
    sphingosine_patterns_full = [
        '[N][C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCC',  # General full backbone with usual stereochemistry
        '[N][C@@H](CO)[C@@H](O)C=CC', # To catch possible other arrangements
    ]
    
    # Check for a generic sphingosine backbone match among possible patterns
    if all(not mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in sphingosine_patterns_full):
        return False, "No sphingosine backbone with appropriate stereochemistry and length found"

    # Define β-D-glucose attachment with flexibility on connecting atom
    glucose_pattern = Chem.MolFromSmarts('OC[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O')
    # Look for patterns where glucose is connected to longer hydrocarbon chains
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No β-D-glucose moiety found attached to the primary hydroxyl"

    # Check for amide bond to fatty acyl chain
    fatty_acyl_amide_pattern = Chem.MolFromSmarts('C(=O)N')
    if not mol.HasSubstructMatch(fatty_acyl_amide_pattern):
        return False, "Amide linkage to fatty acyl chain not found"
    
    # If all patterns match, classify as glucosylceramide
    return True, "Contains sphingosine backbone, amide linkage, and β-D-glucose moiety"