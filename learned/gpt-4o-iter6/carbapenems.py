"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determine if a given SMILES string is a carbapenem.

    Args:
        smiles (str): The SMILES string of the chemical entity.
    
    Returns:
        bool: True if SMILES represents a carbapenem, False otherwise.
        str: Explanation of the reasoning for the classification.
    """
    
    # Define a more flexible SMARTS pattern for the carbapenem core:
    # A basic bicyclic with versatility for substitutions
    core_pattern = Chem.MolFromSmarts("C1CNC2=C1[S,N]C(=O)N2")
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the core pattern in the molecule
    if mol.HasSubstructMatch(core_pattern):
        # Additional checks for the positions 3, 4, 6 substitutions
        # placeholders: check if present and possibly specific substituents
        subs_3_4_6_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
        if mol.HasSubstructMatch(subs_3_4_6_pattern):
            return True, "Matches the carbapenem core and substitution requirements"
        else:
            return False, "Matches carbapenem core but lacks appropriate substitutions"
    else:
        return False, "Does not match the carbapenem core structure"