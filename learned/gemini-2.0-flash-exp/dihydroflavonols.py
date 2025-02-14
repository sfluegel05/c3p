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

    # Define the flavanone core SMARTS with the hydroxyl group at C3.
    # Flexible SMARTS that match correctly.
    flavanone_core_smarts = "[c]1[c]([OH])[C]([c]([c]1)[O])=[O][C]"
    flavanone_core_pattern = Chem.MolFromSmarts(flavanone_core_smarts)
    
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Flavanone core with 3-hydroxyl not found"
    

    return True, "Molecule is a dihydroflavonol"