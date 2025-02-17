"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: enone – alpha,beta-unsaturated ketones (R(1)R(2)C=CR(3)-C(=O)R(4) with R(4) ≠ H)
A molecule containing an enone should contain a conjugated C=C–C(=O)–C fragment.
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    A chemical entity is classified as an enone if it contains an alpha,beta-unsaturated ketone
    motif of general structure R(1)R(2)C=CR(3)-C(=O)R(4) (with R(4) not being H).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is an enone, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the enone substructure.
    # This pattern corresponds to a conjugated C=C–C(=O)–C motif. The carbonyl must be bonded to a carbon (not H)
    enone_pattern = Chem.MolFromSmarts("[CX3]=[CX3]-C(=O)[CX3]")
    if enone_pattern is None:
        return False, "Could not generate SMARTS for enone"

    # Check if the molecule contains the enone substructure
    if mol.HasSubstructMatch(enone_pattern):
        return True, "Contains an alpha,beta-unsaturated ketone (enone) motif: C=C-C(=O)-C"
    else:
        return False, "Does not contain the enone motif (conjugated C=C-C(=O)-C structure)"
        
# Example test runs (uncomment to test)
# examples = [
#     "C=CC(=O)C",  # minimal enone (methyl vinyl ketone derivative)
#     "CC1=CC(=O)CC(C)(C)C1"  # e.g., isophorone
# ]
# for sm in examples:
#     result, reason = is_enone(sm)
#     print(f"SMILES: {sm}\nResult: {result}, Reason: {reason}\n")