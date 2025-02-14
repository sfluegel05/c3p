"""
Classifies: CHEBI:36916 cation
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is a monoatomic or polyatomic species having one or more elementary charges of the proton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cation, False otherwise.
        str: Reason for classification.
    """

    # Handle the simple proton case first
    if smiles == "[H+]":
         return True, "Single proton cation"
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for negative charges first
    for atom in mol.GetAtoms():
         if atom.GetFormalCharge() < 0:
              return False, "Contains negatively charged atoms, not a simple cation"

    has_positive_charge = False
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() > 0:
            # Check for counterions such as Cl, O
            if atom.GetAtomicNum() == 17 or atom.GetAtomicNum() == 8:
                 return False, "Positive charge on inappropriate atom. Likely salt, not cation"

            has_positive_charge = True

    if has_positive_charge:
        # Verify molecule doesn't have other molecules
       
        return True, "Contains at least one positive charge and no negative charges."
    else:
        return False, "No positive charge found."