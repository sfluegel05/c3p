"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: Anthoxanthin – flavonoid pigments in plants that are typically based on a flavone (2-phenylchromen-4-one) scaffold,
although extensive substitution (e.g., methoxylation, hydroxylation, glycosylation) is common.
The new approach uses revised SMARTS that focus on the core connectivity of a flavone while allowing substituents.
"""

from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin (a flavonoid pigment) based on its SMILES string.
    Here we search for a flavone (2-phenylchromen-4-one) core, taking into account substituent tolerance.
    
    We use two SMARTS patterns which represent the core connectivity of the flavone skeleton:
      Pattern 1: "c1ccc(-c2oc(=O)c3ccccc3c2)cc1"
      Pattern 2: "c1ccc(-c2cc(=O)c3occcc3c2)cc1"
    
    These patterns capture the two typical ways in which the 2-phenylchromen-4-one core is drawn.
    If either is found as a substructure, we classify the molecule as an anthoxanthin.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as anthoxanthin, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for the flavone core.
    # Pattern 1: This pattern represents a flavone core where the 2-phenyl group is attached directly
    # to a chromen-4-one system. Substituents on the aromatic rings will not interfere with matching.
    flavone_smarts_1 = "c1ccc(-c2oc(=O)c3ccccc3c2)cc1"
    flavone_pattern_1 = Chem.MolFromSmarts(flavone_smarts_1)
    if flavone_pattern_1 is None:
        return False, "Error in parsing flavone substructure (pattern 1)"
    
    # Pattern 2: An alternate connectivity order for the heterocyclic part.
    flavone_smarts_2 = "c1ccc(-c2cc(=O)c3occcc3c2)cc1"
    flavone_pattern_2 = Chem.MolFromSmarts(flavone_smarts_2)
    if flavone_pattern_2 is None:
        return False, "Error in parsing flavone substructure (pattern 2)"
    
    # Check if either pattern matches the molecule.
    if mol.HasSubstructMatch(flavone_pattern_1) or mol.HasSubstructMatch(flavone_pattern_2):
        return True, "Molecule contains a flavone (2-phenylchromen-4-one) core typical of anthoxanthins."
    else:
        return False, "Molecule does not contain the expected flavone core required for anthoxanthins."

# Example usage (for testing purposes only):
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