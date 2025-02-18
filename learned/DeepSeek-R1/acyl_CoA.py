"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: CHEBI:17996 acyl-CoA
"""
from rdkit import Chem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester formed between coenzyme A's thiol group and a carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the thioester group connected to pantetheine chain (S-C(=O)-C-C-N-C(=O))
    thioester_pantetheine_pattern = Chem.MolFromSmarts("[S](C(=O))C-C-N-C(=O)")
    
    # Define adenine moiety pattern (part of CoA structure)
    adenine_pattern = Chem.MolFromSmarts("[n]1cnc2ncnc12")

    # Check for both substructures
    has_thioester = mol.HasSubstructMatch(thioester_pantetheine_pattern)
    has_adenine = mol.HasSubstructMatch(adenine_pattern)

    if has_thioester and has_adenine:
        return True, "Contains thioester group linked to CoA pantetheine chain and adenine moiety"
    else:
        reasons = []
        if not has_thioester:
            reasons.append("Missing thioester linkage to CoA backbone")
        if not has_adenine:
            reasons.append("Missing adenine moiety characteristic of CoA")
        return False, ", ".join(reasons) if reasons else "Does not match acyl-CoA criteria"