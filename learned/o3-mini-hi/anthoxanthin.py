"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: Anthoxanthin – flavonoid pigments in plants that are typically water-soluble and range in color from colorless/creamy to yellow.
This program checks if the input SMILES string contains a flavone (2-phenylchromen-4-one) core, which is essential for anthoxanthins.
"""

from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin (a type of flavonoid pigment) based on its SMILES string.
    Anthoxanthins are typically flavones (2-phenylchromen-4-one) or their derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is an anthoxanthin (contains a flavone core), False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a basic flavone core:
    # This pattern represents a 2-phenylchromen-4-one scaffold:
    # - "c1ccc(cc1)" identifies a benzene ring (the B ring)
    # - "-c2oc(=O)c3ccccc3o2" identifies the fused heterocyclic part (rings A and C) of the flavone.
    flavone_smarts = "c1ccc(cc1)-c2oc(=O)c3ccccc3o2"
    flavone_pattern = Chem.MolFromSmarts(flavone_smarts)
    if flavone_pattern is None:
        return False, "Error in parsing substructure pattern"
    
    # Check if the molecule contains the flavone core.
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Molecule contains the 2-phenylchromen-4-one (flavone) scaffold typical of anthoxanthins."
    else:
        return False, "Molecule does not contain the expected flavone core required for anthoxanthins."

# Example usage (uncomment to test):
# smiles_examples = [
#     "COc1ccc(cc1OC)-c1cc(=O)c2c(OC)c(OC)c(OC)cc2o1",  # sinensetin
#     "COc1ccc(cc1)-c1oc2c(OC)c(OC)cc(O)c2c(=O)c1O",     # tambulin
#     "COc1c(O)cc(O)c2c1oc(cc2=O)-c1ccccc1",             # wogonin
# ]
#
# for s in smiles_examples:
#     result, reason = is_anthoxanthin(s)
#     print(f"SMILES: {s}\nClassified as anthoxanthin? {result}\nReason: {reason}\n")