"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: CHEBI: 3α-hydroxy steroid (heuristic)
A 3α-hydroxy steroid is defined as a steroid having a typical fused tetracyclic (three six-membered rings and one five-membered ring) core,
with a hydroxyl (-OH) substituent on one of the six‐membered rings of that core.
This heuristic uses a generalized SMARTS for the steroid nucleus and looks for a [C;R6][OX2H] pattern.
Note: We relax the stereochemistry matching (useChirality=False) so that variations in SMILES representation do not miss a valid steroid core.
"""

from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3α-hydroxy steroid based on its SMILES string using heuristic substructure matches.
    
    The following steps are taken:
      1. A more general steroid nucleus is defined as a fused tetracyclic ring system (three six-membered rings and one five-membered ring)
         using the SMARTS "C1CCC2C3CCC4C(C3CCC2C1)CCC4". We relax stereochemistry by using useChirality=False.
      2. A hydroxyl group on a six‐membered ring carbon is searched for using the SMARTS "[C;R6][OX2H]".
         At least one such match must fall on an atom that is part of the steroid nucleus.
    
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
    
    # Define a more general SMARTS for a steroid nucleus (fused tetracyclic core)
    steroid_core_smarts = "C1CCC2C3CCC4C(C3CCC2C1)CCC4"
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
        return False, "Error generating steroid core pattern"
    
    # Use relaxed stereo matching (useChirality=False) so that variations in stereochemistry do not prevent a match.
    if not mol.HasSubstructMatch(steroid_core, useChirality=False):
        return False, "Steroid nucleus not found based on the heuristic core pattern"
    
    # Retrieve all atom indices that are the part of a steroid core match.
    core_matches = mol.GetSubstructMatches(steroid_core, useChirality=False)
    core_atoms = set()
    for match in core_matches:
        core_atoms.update(match)
    if not core_atoms:
        return False, "Steroid core atoms could not be identified"
    
    # Define a SMARTS pattern for a hydroxyl group (-OH) on a six-membered ring carbon.
    oh_pattern = Chem.MolFromSmarts("[C;R6][OX2H]")
    if oh_pattern is None:
        return False, "Error generating hydroxyl group pattern"
    
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if not oh_matches:
        return False, "No hydroxyl group found on a six-membered ring carbon in the molecule"
    
    # Verify that at least one hydroxyl-bearing carbon is part of the steroid nucleus.
    for match in oh_matches:
        c_atom_idx = match[0]
        if c_atom_idx in core_atoms:
            return True, "Molecule contains a steroid nucleus with a hydroxyl group on a six-membered ring consistent with a 3α-hydroxy steroid"
    
    return False, "Hydroxyl group on a steroid nucleus not found"

# Example usage:
if __name__ == "__main__":
    # Example test SMILES from one of the provided examples:
    test_smiles = "C1[C@@]2([C@@]([C@@]3([C@](C[C@H](O)CC3)(C1)[H])C)(CC[C@@]4([C@H](CC[C@@]24[H])OS(O)(=O)=O)C)[H])[H]"
    result, reason = is_3alpha_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)