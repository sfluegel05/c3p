"""
Classifies: CHEBI:16158 steroid sulfate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a sulfuric ester obtained by the condensation of a hydroxy group of any steroid with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid backbone pattern using SMARTS - a tetracyclic carbon ring structure
    steroid_pattern = Chem.MolFromSmarts("C12CCC3C(C1)CCC4C(C3)CCC4")
    
    # Check if the steroid backbone exists in the molecule
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Define sulfate group pattern using SMARTS
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O-]")

    # Find all sulfate matches
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    if not sulfate_matches:
        return False, "No sulfate groups found"

    # Check for connectivity: sulfate group (ester) must be connected to the steroid at hydroxy site
    # Define generic ester group pattern
    ester_sulfate_pattern = Chem.MolFromSmarts("C-O-S(=O)(=O)[O-]")
    
    # Search for the hydroxy group on steroid modified by sulfate
    for match in sulfate_matches:
        # Extract the atom indices of the match (we want to ensure connection to steroid backbone)
        if any(mol.HasSubstructMatch(ester_sulfate_pattern)):
            return True, "Contains sulfate ester linked to steroid backbone at hydroxy group"
    
    return False, "Sulfate groups found, but not esterified to steroid backbone"

# Test with an example structure
smiles = "C1=C2C(CC[C@]3([C@@]4(CC[C@@H]([C@]4(CC[C@@]32[H])C)O)[H])[H])=CC(=C1)OS(O)(=O)=O"
result, reason = is_steroid_sulfate(smiles)
print(f"Result: {result}, Reason: {reason}")