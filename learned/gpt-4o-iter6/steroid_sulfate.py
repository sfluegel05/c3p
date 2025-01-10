"""
Classifies: CHEBI:16158 steroid sulfate
"""
from rdkit import Chem

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

    # Define a rough steroid backbone pattern (acknowledge complexity and why exact match is challenging)
    # The SMARTS string "C1CCC2C(C1)CCC3C2CCC4C3CCC4" represents a basic cholesterol steroid ring pattern, not too specific for all steroids due to variability.

    # Check for presence of a rough steroid backbone
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CCC4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Look for sulfate groups esters
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)O")  # Corrected sulfate ester bond

    # Check if the sulfate group exists and is connected to the steroid backbone
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    if not sulfate_matches:
        return False, "No sulfate groups connected as esters"

    # Check connectivity to alcohols (which would have been the type to esterify)
    hydroxyl_pattern = Chem.MolFromSmarts("CO")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups found for esterification"

    for match in sulfate_matches:
        # Ensure ester linkage to the hydroxy site on predicted steroid backbone
        atom_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'O']
        if any(mol.GetAtomWithIdx(idx).GetNeighbors()[0].GetSymbol() == 'C' for idx in atom_indices):
            return True, "Contains sulfate ester linked to steroid backbone at a hydroxy group"

    return False, "Sulfate groups found but not linked properly to steroid backbone"

# Example structure for testing (estrone sodium sulfate)
smiles = "[Na].C[C@]12CC[C@H]3C(=CCc4cc(OS([O-])(=O)=O)ccc34)[C@@H]1CCC2=O"
result, reason = is_steroid_sulfate(smiles)
print(f"Result: {result}, Reason: {reason}")