"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: CHEBI:77834 3alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid is a steroid with a hydroxy group at the 3-position in the alpha configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C@]1(CCC[C@@H]2[C@@]3([H])CCC[C@]3([H])C[C@@H]2[C@@]1([H])C)")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for 3alpha-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@@]([H])(O)([H])")
    if not any(mol.GetAtomWithIdx(idx).HasProp("_3alphaHydroxyPos") for idx in mol.GetSubstructMatches(hydroxy_pattern)):
        return False, "No 3alpha-hydroxy group found"
    
    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 600:
        return False, "Molecular weight outside typical range for steroids"
    
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 3 or n_rings > 5:
        return False, "Number of rings outside typical range for steroids"
    
    return True, "Contains a steroid backbone with a hydroxy group at the 3-position in the alpha configuration"