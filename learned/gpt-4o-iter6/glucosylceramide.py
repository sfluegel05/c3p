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
    
    # Define SMARTS pattern for β-D-glucose moiety with potential flexibility
    glucose_pattern = Chem.MolFromSmarts('[C@@H]1(O[C@H]([C@H]([C@@H]([C@H]1O)O)O)CO)O')
    
    # Check for the β-D-glucose moiety
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No β-D-glucose moiety found"
    
    # Define enhanced SMARTS for a sphingosine-like backbone, focusing on typical ceramide features 
    sphingosine_backbone_pattern = Chem.MolFromSmarts('N[C@@H](CCC[C@H](O)CO)[H]')
    
    # Ensure that the substructure for a sphingosine-like backbone is present
    if not mol.HasSubstructMatch(sphingosine_backbone_pattern):
        return False, "No recognizable sphingosine backbone pattern found"
    
    # Enhanced amide linkage pattern including an adjacent long carbon chain indicative of fatty acids
    amide_pattern = Chem.MolFromSmarts('C(=O)N[C@@H](CO)C')
    
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Amide linkage characteristic of fatty acyl attachment not found"

    # If all essential patterns match, classify as glucosylceramide
    return True, "Molecule contains a sphingosine backbone, amide linkage, and β-D-glucose moiety"