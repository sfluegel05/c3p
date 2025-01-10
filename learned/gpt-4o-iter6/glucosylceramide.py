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

    # Redefine SMARTS patterns with accurate stereochemistry for β-D-glucose moiety
    glucose_pattern = Chem.MolFromSmarts('[C@@H]1O[C@H]([C@@H](O)[C@H](O)[C@H]([C@H]1O)O)CO')

    # Check for the glucosyl moiety
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No β-D-glucose moiety found attached to the primary hydroxyl"

    # Flexible sphingosine backbone pattern capturing typical arrangements
    sphingosine_patterns = [
        'N[C@H](CO)C(O)[C@H](O)CCCCCCCCCC', # Typical structure
        'NC[CH2]C(O)[C@H]([CH2]C=C)C' # Possible variations in connection
    ]

    if not any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in sphingosine_patterns):
        return False, "No compatible sphingosine backbone pattern found"

    # Amide bond detection for fatty acyl chain
    fatty_acyl_amide_pattern = Chem.MolFromSmarts('C(=O)N')
    if not mol.HasSubstructMatch(fatty_acyl_amide_pattern):
        return False, "Amide linkage to fatty acyl chain not found"
    
    # If all essential patterns match, classify as glucosylceramide
    return True, "Molecule contains a sphingosine backbone, amide linkage, and β-D-glucose moiety"