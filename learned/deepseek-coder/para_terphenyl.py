"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: CHEBI:XXXXX para-terphenyl
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl based on its SMILES string.
    A para-terphenyl is a molecule with a 1,4-diphenylbenzene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the 1,4-diphenylbenzene pattern
    para_terphenyl_pattern = Chem.MolFromSmarts("c1ccc(-c2ccccc2)cc1-c3ccccc3")
    
    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(para_terphenyl_pattern):
        return True, "Contains a 1,4-diphenylbenzene skeleton"
    else:
        return False, "Does not contain a 1,4-diphenylbenzene skeleton"