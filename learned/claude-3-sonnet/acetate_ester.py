"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: CHEBI:35914 acetate ester
An acetate ester is any carboxylic ester where the carboxylic acid component is acetic acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Look for ester pattern (-O-C(=O)-) with flexible SMARTS
    ester_pattern = Chem.MolFromSmarts("[OX2]C(=O)[A;!$(OC)]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check if any ester group is an acetate ester
    for ester_match in ester_matches:
        ester_atoms = [mol.GetAtomWithIdx(idx) for idx in ester_match]
        carbonyl_atom = ester_atoms[1]
        
        # Look for acetate group (-O-C(=O)C) connected to the ester carbonyl
        acetate_pattern = Chem.MolFromSmarts("[OX2]C(=O)(C)[A]")
        acetate_match = mol.GetSubstructMatches(acetate_pattern, atomsToUse=ester_atoms)
        
        if acetate_match:
            return True, "Contains an acetate ester group (-O-C(=O)C)"

    return False, "No acetate ester group found"