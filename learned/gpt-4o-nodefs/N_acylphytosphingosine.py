"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
from rdkit import Chem

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine typically has a phytosphingosine backbone with an acyl group
    attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the N-acylphytosphingosine backbone pattern:
    # − Include hydroxylated long chain (typically 18 carbons)
    # − Carbon chain with secondary amine
    # − Multiple hydroxyl groups, generally two close to the amine group

    # Define SMARTS pattern for phytosphingosine backbone
    # General structure: C(CO)O-CH-(CH2)n-CH(NC)-CH2O
    backbone_pattern = Chem.MolFromSmarts("[C@@H](O)[C@@H](N)[C@@H](O)")  # Looking for amino-alcohol pattern
    
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No phytosphingosine backbone found"

    # Define pattern for an acyl group attached to the nitrogen
    acyl_pattern = Chem.MolFromSmarts("NC(=O)C")  # Acyl group pattern; nitrogen connected to carbonyl
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group attached to the nitrogen"

    return True, "Contains phytosphingosine backbone with N-acyl group"

# Test the function
example_smiles = "CCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC"
result, reason = is_N_acylphytosphingosine(example_smiles)
print(f"Result: {result}, Reason: {reason}")