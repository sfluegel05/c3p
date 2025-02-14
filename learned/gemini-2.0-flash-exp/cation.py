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

    #split into molecules, handle each separately
    smiles_list = smiles.split(".")
    
    is_cation = False
    reason = "Not a cation"

    for s in smiles_list:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
           return False, "Invalid SMILES string"

        has_positive_charge = False
        for atom in mol.GetAtoms():
            if atom.GetFormalCharge() > 0:
                has_positive_charge = True
                break # No need to check further atoms in this molecule

        if has_positive_charge:
          is_cation = True
          reason = "Contains at least one positive charge"

    return is_cation, reason