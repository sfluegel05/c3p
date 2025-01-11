"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: CHEBI:29278 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound with a thiol group (-SH) attached to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the thiol group pattern (-SH attached to a carbon, aliphatic or aromatic)
    thiol_pattern = Chem.MolFromSmarts("[C,c][SH1]")
    
    # Check if the molecule contains the thiol pattern
    if mol.HasSubstructMatch(thiol_pattern):
        # Exclude molecules where sulfur is part of other functional groups
        sulfide_pattern = Chem.MolFromSmarts("[C,c][S][C,c]")
        disulfide_pattern = Chem.MolFromSmarts("[C,c][S][S][C,c]")
        sulfoxide_pattern = Chem.MolFromSmarts("[C,c][S](=O)[C,c]")
        sulfone_pattern = Chem.MolFromSmarts("[C,c][S](=O)(=O)[C,c]")
        thiocarbonyl_pattern = Chem.MolFromSmarts("[C,c][S]=[C,c]")
        
        if (not mol.HasSubstructMatch(sulfide_pattern) and
            not mol.HasSubstructMatch(disulfide_pattern) and
            not mol.HasSubstructMatch(sulfoxide_pattern) and
            not mol.HasSubstructMatch(sulfone_pattern) and
            not mol.HasSubstructMatch(thiocarbonyl_pattern)):
            return True, "Contains a thiol group (-SH) attached to a carbon atom"
        else:
            return False, "Sulfur is part of a non-thiol functional group"
    else:
        return False, "No thiol group (-SH) found attached to a carbon atom"