"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester contains the substructure CH3C(=O)O-.

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

    # Define SMARTS pattern for acetate group (CH3C(=O)O-)
    acetate_pattern = Chem.MolFromSmarts("CC(=O)O")

    # Search for the substructure
    matches = mol.GetSubstructMatches(acetate_pattern)

    # Check if the number of matches is exactly 1
    if len(matches) == 1:
       # Check if the molecule is large enough
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 100:
             return False, "Molecule too small to be an acetate ester."
        return True, "Contains an acetate group"
    elif len(matches) > 1:
        return False, f"Contains {len(matches)} acetate groups, needs exactly 1"
    else:
        return False, "Does not contain an acetate group"