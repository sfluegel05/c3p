"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: CHEBI: 3α-hydroxy steroid (heuristic)
A 3α-hydroxy steroid is defined as a steroid (i.e. containing a typical fused tetracyclic system of three six-membered rings and one five-membered ring)
with a hydroxyl (-OH) substituent on a ring carbon in what is presumed to be the A‐ring. 
Since automatic numbering and precise stereochemistry assignment is challenging, the approach is to (1) detect a steroid nucleus by a common core SMARTS,
and then (2) require at least one –OH on a six‐membered ring atom within that fused system. 

This heuristic may fail on borderline cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3α-hydroxy steroid based on its SMILES string using heuristic SMARTS substructure
    matches:
      1. Checks that the molecule contains a steroid core (a fused tetracyclic system resembling a cyclopentanoperhydrophenanthrene).
      2. Checks that at least one six-membered ring carbon in that core carries an –OH substituent.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3α-hydroxy steroid, False otherwise.
        str: A message with the reason for the classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Heuristic steroid nucleus SMARTS.
    # This pattern looks for a fused ring system reminiscent of the steroid core. 
    # (Note: many steroid cores exist; this is one commonly used pattern.)
    steroid_core_smarts = "C1CC2CCC3C(C2)CC[C@]13"
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Steroid nucleus not found based on the heuristic core pattern"
    
    # Get all matching atom indices for the steroid core pattern.
    # (There may be more than one match; we take the union of atoms from all matches.)
    core_matches = mol.GetSubstructMatches(steroid_core)
    core_atoms = set()
    for match in core_matches:
        core_atoms.update(match)
    if len(core_atoms) == 0:
        return False, "Steroid core atoms could not be identified"
    
    # Define a heuristic SMARTS for an aliphatic hydroxyl group attached to a ring carbon
    # in a six-membered ring. The "[C;R6]" restricts the carbon to be in a six-membered ring.
    oh_pattern = Chem.MolFromSmarts("[C;R6][OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if len(oh_matches) == 0:
        return False, "No hydroxyl group found on a six-membered ring carbon"
    
    # Check that at least one hydroxyl-bearing carbon is part of the steroid nucleus.
    found = False
    for match in oh_matches:
        c_idx = match[0]  # The carbon atom in the [C;R6][OH] pattern.
        if c_idx in core_atoms:
            found = True
            break
            
    if not found:
        return False, "Hydroxyl group not found on the identified steroid core"
    
    # (Optional additional criteria can be added here if further refinement is desired,
    #     such as molecular weight or the number of rings, etc.)
    
    return True, "Molecule contains a steroid nucleus with a hydroxyl group on a six-membered ring consistent with a 3α-hydroxy steroid"

# Example usage (you can remove or comment out these lines in production):
if __name__ == "__main__":
    # One example SMILES from the provided list.
    example_smiles = "C1[C@@]2([C@@]([C@@]3([C@](C[C@H](O)CC3)(C1)[H])C)(CC[C@@]4([C@H](CC[C@@]24[H])OS(O)(=O)=O)C)[H])[H]"
    result, reason = is_3alpha_hydroxy_steroid(example_smiles)
    print("Result:", result)
    print("Reason:", reason)