"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: Anthoxanthin – flavonoid pigments in plants that range in color from colorless/creamy to yellow.
Anthoxanthins are generally based on a flavone (2-phenylchromen-4-one) scaffold, although substitutions
may lead to variations in the drawn structure.
The program uses two SMARTS patterns to capture the typical flavone core connectivity.
"""

from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin (a flavonoid pigment) based on its SMILES string.
    We search for a flavone (2-phenylchromen-4-one) core, but because of connectivity variations due
    to substituents, two patterns are used.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an anthoxanthin (contains a flavone core), False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for the flavone core.
    # Pattern 1: assumes the oxygen in the heterocycle comes immediately after the bond
    # from the B ring. This pattern has been used previously:
    flavone_smarts_1 = "c1ccc(cc1)-c2oc(=O)c3ccccc3o2"
    flavone_pattern_1 = Chem.MolFromSmarts(flavone_smarts_1)
    if flavone_pattern_1 is None:
        return False, "Error in parsing flavone substructure (pattern 1)"
    
    # Pattern 2: allows for an alternate ordering in the heterocycle,
    # where the carbonyl and the oxygen are interchanged.
    flavone_smarts_2 = "c1ccc(cc1)-c2cc(=O)c3occcc23"
    flavone_pattern_2 = Chem.MolFromSmarts(flavone_smarts_2)
    if flavone_pattern_2 is None:
        return False, "Error in parsing flavone substructure (pattern 2)"

    # Check if either pattern matches the molecule.
    match1 = mol.HasSubstructMatch(flavone_pattern_1)
    match2 = mol.HasSubstructMatch(flavone_pattern_2)
    
    if match1 or match2:
        return True, "Molecule contains a flavone (2-phenylchromen-4-one) core typical of anthoxanthins."
    else:
        return False, "Molecule does not contain the expected flavone core required for anthoxanthins."

# Example usage:
if __name__ == '__main__':
    smiles_examples = [
        "COc1ccc(cc1OC)-c1cc(=O)c2c(OC)c(OC)c(OC)cc2o1",  # sinensetin
        "COc1ccc(cc1)-c1oc2c(OC)c(OC)cc(O)c2c(=O)c1O",     # tambulin
        "COc1c(O)cc(O)c2c1oc(cc2=O)-c1ccccc1",             # wogonin
        "COc1c(O)c2c(oc(cc2=O)-c2ccc(O)cc2)c(OC)c1OC",      # xanthomicrol
    ]
    
    for smi in smiles_examples:
        result, reason = is_anthoxanthin(smi)
        print(f"SMILES: {smi}\nClassified as anthoxanthin? {result}\nReason: {reason}\n")