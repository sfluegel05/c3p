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
    
    # Broad pattern for phytosphingosine backbone recognition:
    # Allow for multiple stereochemistry configurations and variations
    # Include hydroxylated long chain and secondary amide connections
    # Pattern components: Long chain with hydroxyls and a secondary amine
    backbone_pattern = Chem.MolFromSmarts("[C@@H,O][C@@H](N(C=O))([C@@H,O])+") # A rougher generalization
    
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No N-acylphytosphingosine backbone pattern found"

    # Acyl group should be attached to the nitrogen, we check for a carbonyl group connected to the nitrogen
    acyl_group_pattern = Chem.MolFromSmarts("N[C;$(C=O)]")  # Highlight carbonyl linkage
    
    if not mol.HasSubstructMatch(acyl_group_pattern):
        return False, "No acyl group attached to the nitrogen"

    return True, "Contains N-acylphytosphingosine pattern with a matching backbone and acyl group"