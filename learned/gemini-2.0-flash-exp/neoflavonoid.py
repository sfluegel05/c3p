"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is a 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the chromene core with a substituent at C4 position
    # This SMARTS explicitly defines the double bonds and the carbons in the core structure, and the attachment point for the aryl substituent.
    core_pattern = Chem.MolFromSmarts('c1ccccc1-c2-O-C(=[C])-[C](-[*])-C=C2')
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Does not contain a 1-benzopyran core."
    
    # Define the aryl substituent as a generic substituted or unsubstituted phenyl ring, connected via a single bond to the core
    #   and the SMARTS includes the connection to the carbon of the core.
    aryl_substituent_pattern = Chem.MolFromSmarts('c1ccccc1-[*]')
    if not mol.HasSubstructMatch(aryl_substituent_pattern):
        return False, "Does not have an aryl substituent at position 4 of the 1-benzopyran core."
    
    # Combined pattern to check for the core and the aryl substituent attached at position 4
    combined_pattern = Chem.MolFromSmarts('c1ccccc1-c2-O-C(=[C])-[C](-c3ccccc3)-C=C2')
    if not mol.HasSubstructMatch(combined_pattern):
         return False, "Does not contain a 1-benzopyran core with an aryl substituent at position 4."


    return True, "Contains a 1-benzopyran core with an aryl substituent at position 4."