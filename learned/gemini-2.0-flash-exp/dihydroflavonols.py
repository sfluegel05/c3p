"""
Classifies: CHEBI:48039 dihydroflavonols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is a flavanone with a hydroxyl group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavanone core SMARTS with the hydroxyl group at C3
    # Numbering the flavanone system so that we can uniquely identify the 3-hydroxyl.
    # The key is to have a C=O at position 4 and an -OH at position 3 of the fused rings
    flavanone_core_smarts = "[C]1[c]2[c]([c]([OH])[C](=[O])[C]1)[cH][cH][cH][cH]2"
    flavanone_core_pattern = Chem.MolFromSmarts(flavanone_core_smarts)
    
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Flavanone core with 3-hydroxyl not found"
    

    return True, "Molecule is a dihydroflavonol"