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

    # Define SMARTS pattern for β-D-glucose moiety with proper stereochemistry
    glucose_pattern = Chem.MolFromSmarts('OC[C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O')
    
    # Check for the glucosyl moiety
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No β-D-glucose moiety found attached to the primary hydroxyl"

    # Refined sphingosine backbone pattern that captures typical arrangements
    sphingosine_pattern = Chem.MolFromSmarts('N[C@@H](CO[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](O)1)C')
    
    # Ensure that the substructure for a sphingosine-like backbone is present
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No compatible sphingosine backbone pattern found"

    # Amide bond detection for linking the fatty acyl chain (simple representation)
    fatty_acyl_amide_pattern = Chem.MolFromSmarts('C(=O)N')
    
    if not mol.HasSubstructMatch(fatty_acyl_amide_pattern):
        return False, "Amide linkage to fatty acyl chain not found"
    
    # If all essential patterns match, classify as glucosylceramide
    return True, "Molecule contains a sphingosine backbone, amide linkage, and β-D-glucose moiety"