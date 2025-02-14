"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: CHEBI:35523 acetate ester
An acetate ester is any carboxylic ester where the carboxylic acid component is acetic acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for acetate ester substructure pattern
    acetate_ester_pattern = Chem.MolFromSmarts("CC(=O)OC")
    acetate_ester_matches = mol.GetSubstructMatches(acetate_ester_pattern)

    # Ensure there is only one acetate ester group
    if len(acetate_ester_matches) == 1:
        # Check for additional structural requirements
        # (e.g., no other ester bonds, specific molecular weight range, etc.)
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if 50 < mol_wt < 500:  # Typical molecular weight range for acetate esters
            return True, "Contains an acetate group as part of an ester bond, and meets additional structural requirements"

    return False, "Does not contain an acetate ester group, or does not meet additional structural requirements"