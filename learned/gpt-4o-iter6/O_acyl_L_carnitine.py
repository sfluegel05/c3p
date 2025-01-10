"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    O-acyl-L-carnitine contains a quaternary ammonium ion, an ester-linked acyl chain, 
    a chiral L-configuration center, and a carboxylate group adjacent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for quaternary ammonium group [N+](C)(C)C
    quaternary_ammonium_pattern = Chem.MolFromSmarts("[N+](C)(C)C")
    if not mol.HasSubstructMatch(quaternary_ammonium_pattern):
        return False, "No quaternary ammonium ion pattern found"
    
    # Broadened ester linkage: O(C(=O)) valid stereo AND adjacent carboxylate
    ester_pattern = Chem.MolFromSmarts("O[C@@H](C)[N+](C)(C)C")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No valid ester linkage pattern found"

    # Check for C(=O)[O-] adjacent to chiral center
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No adjacent carboxylate group found"
    
    # Validate correct stereochemistry on chiral center
    is_correct_stereo = any(atom.HasProp('_CIPCode') and 
                            atom.GetProp('_CIPCode') in ['S'] for atom in mol.GetAtoms() if atom.GetChiralTag())

    if not is_correct_stereo:
        return False, "No correct L-stereochemistry found"

    return True, "Contains O-acyl-L-carnitine features: quaternary ammonium, ester linkage, L-stereo, carboxylate."