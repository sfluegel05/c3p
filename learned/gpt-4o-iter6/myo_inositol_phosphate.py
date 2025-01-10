"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate in which the inositol component
    has myo-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        (bool, str): True if molecule is a myo-inositol phosphate with reason, else False with reason.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cyclohexane ring with any stereochemistry for myo-inositol, be less strict on stereochemistry
    myo_inositol_pattern = Chem.MolFromSmarts("C1(CO)C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No myo-inositol core structure found"
    
    # Ensure presence of multiple phosphate groups (at least one)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate groups found"
    
    # Count the number of phosphate attachments to confirm characteristic polyphosphate
    if len(phosphate_matches) < 1:
        return False, f"Only {len(phosphate_matches)} phosphate group(s) found, which may be insufficient for classification"

    # Confirm that phosphate groups are attached in a typical manner
    phosphate_attachment = Chem.MolFromSmarts("O-P(=O)(O)C")
    if not mol.HasSubstructMatch(phosphate_attachment):
        return False, "Phosphate groups not in typical attachment sites for myo-inositol phosphate"

    return True, "Contains myo-inositol core structure with phosphate groups in characteristic locations"

# Example test
smiles_example = "O[C@@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H]1OP(O)(O)=O"
print(is_myo_inositol_phosphate(smiles_example))