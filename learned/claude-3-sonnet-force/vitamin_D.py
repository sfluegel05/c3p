"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: CHEBI:27279 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D based on its SMILES string.
    Vitamin D molecules are fat-soluble hydroxy seco-steroids with a cholesterol-like backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cholesterol-like backbone (4 fused rings, alkyl side chain)
    steroid_pattern = Chem.MolFromSmarts("[*]1[*]2[*]3[*]4[*]1[*]2[*]3[*]4[CH2]CC")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Look for seco-structure (broken ring)
    seco_pattern = Chem.MolFromSmarts("[*]1[*]2[*]3[*]4[*]1[*]2[*]3=[*]4")
    if not mol.HasSubstructMatch(seco_pattern):
        return False, "No seco-structure found"

    # Look for hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) < 1:
        return False, "No hydroxy group found"

    # Count rotatable bonds to verify alkyl chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Alkyl chain too short for vitamin D"

    # Check molecular weight - vitamin D typically 300-500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 500:
        return False, "Molecular weight out of range for vitamin D"

    return True, "Contains cholesterol-like steroid backbone with seco-structure and hydroxy group"