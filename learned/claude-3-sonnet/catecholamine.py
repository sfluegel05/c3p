"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: CHEBI:38469 catecholamine
4-(2-Aminoethyl)pyrocatechol [4-(2-aminoethyl)benzene-1,2-diol] and derivatives formed by substitution.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        str: Classification of the molecule as a catecholamine
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES string"

    # Core catecholamine structure: pyrocatechol ring with 2-aminoethyl side chain
    core_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4][CX4]c1c(O)c(O)cc(c1)O")
    if mol.HasSubstructMatch(core_pattern):
        
        # Check for specific catecholamines
        dopamine_pattern = Chem.MolFromSmarts("NCCc1ccc(O)c(O)c1")
        if mol.HasSubstructMatch(dopamine_pattern):
            return "Dopamine"
        
        noradrenaline_pattern = Chem.MolFromSmarts("[C@@H](N)(O)Cc1ccc(O)c(O)c1")
        if mol.HasSubstructMatch(noradrenaline_pattern):
            return "(R)-Noradrenaline"
        noradrenaline_pattern = Chem.MolFromSmarts("[C@H](N)(O)Cc1ccc(O)c(O)c1")
        if mol.HasSubstructMatch(noradrenaline_pattern):
            return "(S)-Noradrenaline"
        
        adrenaline_pattern = Chem.MolFromSmarts("CNC[C@@H](O)c1ccc(O)c(O)c1")
        if mol.HasSubstructMatch(adrenaline_pattern):
            return "(R)-Adrenaline"
        adrenaline_pattern = Chem.MolFromSmarts("CNC[C@H](O)c1ccc(O)c(O)c1")
        if mol.HasSubstructMatch(adrenaline_pattern):
            return "(S)-Adrenaline"
        
        # Check for substituted or complex catecholamines
        sub_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4][CX4]c1c(O)c(O)cc(c1)[X1]")
        if mol.HasSubstructMatch(sub_pattern):
            return "Substituted catecholamine"
        
        return "Primary catecholamine"
    
    return "Not a catecholamine"