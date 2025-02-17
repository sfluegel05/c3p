"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: Anthocyanidin cation
Definition:
  Any organic cation that is an aglycon of anthocyanin cation; they are oxygenated 
  derivatives of flavylium (2-phenylchromenylium).

This script provides the function is_anthocyanidin_cation that takes a SMILES string as input 
and returns a boolean value along with a reason for the classification.
"""
from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    The criteria include:
      - The molecule must be valid.
      - It must be a cation (contain at least one atom with positive formal charge).
      - It should possess a polycyclic ring system (at least three rings).
      - It must contain a core flavylium skeleton (a simplified 2-phenylbenzopyrylium substructure
        with a positively charged oxygen).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an anthocyanidin cation, False otherwise.
        str: Reason for the classification.
    """
    # Attempt to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a positive formal charge (organic cation)
    positive_found = any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms())
    if not positive_found:
        return False, "Molecule is not a cation (no positive formal charge found)"
    
    # Check if the molecule has at least three rings (expected for a flavonoid skeleton)
    if Chem.GetSSSR(mol) < 3:
        return False, "Molecule does not contain enough rings for an anthocyanidin skeleton"
    
    # Define a simplified SMARTS pattern for the flavylium core.
    # This pattern attempts to capture a 2-phenylbenzopyrylium substructure:
    #   - A benzopyrylium ring system where the oxygen is positively charged ([O+])
    #   - A phenyl substituent attached at the 2-position.
    flavylium_smarts = "c1cc2[O+](ccc2)c(c1)-c3ccccc3"
    flavylium = Chem.MolFromSmarts(flavylium_smarts)
    if not mol.HasSubstructMatch(flavylium):
        return False, "Flavylium core pattern not found"
    
    # If all checks pass, we classify the molecule as an anthocyanidin cation.
    return True, "Molecule matches anthocyanidin cation criteria (organic cation with flavylium core)"

# Example usage:
if __name__ == "__main__":
    # A few example SMILES strings (the ones in the problem statement represent glycosides).
    # They may include additional sugar fragments, but our function checks for the underlying
    # flavylium core and the positive charge.
    test_smiles = [
        "COc1cc(cc(OC)c1O)-c1[o+]c2cc(O)cc(O)c2cc1O",  # simplified structure for malvidin
        "OC[C@H]1O[C@@H](Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)c(O)c2)[C@H](O)[C@@H](O)[C@@H]1O"  # one of the examples
    ]
    
    for s in test_smiles:
        result, reason = is_anthocyanidin_cation(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")