"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: CHEBI: 3α-hydroxy steroid (heuristic)
A 3α-hydroxy steroid is defined as a steroid (i.e. having a typical fused tetracyclic system of three six-membered rings and one five-membered ring)
with a hydroxyl (-OH) substituent on a six‐membered ring carbon within the fused system.
Due to variations in stereochemistry notation and extra substituents, this heuristic relies on a generalized SMARTS for the steroid core
and a SMARTS for an aliphatic hydroxyl on a six-membered ring.
"""

from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3α-hydroxy steroid based on its SMILES string using heuristic substructure matches:
      1. It checks that the molecule contains a steroid nucleus – a fused tetracyclic ring system (three six-membered rings and one five-membered ring).
         To avoid strict stereochemical requirements we use a generalized SMARTS pattern.
      2. It then checks that at least one six-membered ring carbon in that core carries an –OH group.
         The pattern [C;R6][OX2H] is used to capture a hydroxyl-bearing carbon that is in a six-membered ring.
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3α-hydroxy steroid, False otherwise.
        str: A message explaining the reasoning behind the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Heuristic SMARTS for a generalized steroid core.
    # This pattern is designed to match a fused tetracyclic system (three six-membered rings and one five-membered ring) 
    # without enforcing specific stereochemistry.
    steroid_core_smarts = "C1CC2CCC3CC4CCCC1C4C23"
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
        return False, "Error in generating steroid core pattern"
    
    # Check if the molecule has a match for the steroid core pattern.
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Steroid nucleus not found based on the heuristic core pattern"
    
    # Get all atom indices that are part of the steroid core match(es)
    core_matches = mol.GetSubstructMatches(steroid_core)
    core_atoms = set()
    for match in core_matches:
        core_atoms.update(match)
    if not core_atoms:
        return False, "Steroid core atoms could not be identified"
    
    # Define a SMARTS for a hydroxyl group on a six-membered ring carbon.
    # [C;R6] means a carbon in a six-membered ring and [OX2H] is a hydroxyl group.
    oh_pattern = Chem.MolFromSmarts("[C;R6][OX2H]")
    if oh_pattern is None:
        return False, "Error in generating hydroxyl group pattern"
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if not oh_matches:
        return False, "No hydroxyl group found on a six-membered ring carbon in the molecule"
    
    # Verify that at least one hydroxyl-bearing carbon is part of the identified steroid nucleus.
    for match in oh_matches:
        c_atom_idx = match[0]  # the carbon atom in the [C;R6][OX2H] match
        if c_atom_idx in core_atoms:
            return True, "Molecule contains a steroid nucleus with a hydroxyl group on a six-membered ring consistent with a 3α-hydroxy steroid"
    
    return False, "Hydroxyl group not found on the identified steroid core"

# Example usage:
if __name__ == "__main__":
    # Test with one example SMILES from the provided list:
    test_smiles = "C1[C@@]2([C@@]([C@@]3([C@](C[C@H](O)CC3)(C1)[H])C)(CC[C@@]4([C@H](CC[C@@]24[H])OS(O)(=O)=O)C)[H])[H]"
    result, reason = is_3alpha_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)