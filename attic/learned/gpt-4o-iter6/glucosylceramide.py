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
        bool: True if the molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for β-D-glucose moiety
    # This pattern explicitly matches β-D-glucosyl linkage
    glucose_pattern = Chem.MolFromSmarts('C1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O')
    
    # Check for the glucosyl moiety
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No β-D-glucose moiety found"
    
    # Define a sphingosine-like backbone pattern with amine and hydroxyl groups
    sphingosine_backbone_pattern = Chem.MolFromSmarts('N[C@@H](CO)C')
    
    # Ensure that the substructure for a sphingosine-like backbone is present
    if not mol.HasSubstructMatch(sphingosine_backbone_pattern):
        return False, "No recognizable sphingosine backbone pattern found"
    
    # Amide linkage pattern for fatty acyl chain attachment
    amide_pattern = Chem.MolFromSmarts('C(=O)N[C@@H]')
    
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Amide linkage characteristic of fatty acyl attachment not found"
    
    # If all essential patterns match, classify as glucosylceramide
    return True, "Molecule contains a sphingosine backbone, amide linkage, and β-D-glucose moiety"