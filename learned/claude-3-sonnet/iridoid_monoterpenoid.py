"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: CHEBI:27894 iridoid monoterpenoid

One of a class of monoterpenoids biosynthesized from isoprene and often intermediates in the biosynthesis of alkaloids. 
Iridoids usually consist of a cyclopentane ring fused to a six-membered oxygen heterocycle; cleavage of a bond in the cyclopentane ring gives rise to the subclass known as secoiridoids.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for core iridoid scaffold
    iridoid_pattern = Chem.MolFromSmarts("[C@H]1[C@H]2[C@H]([C@H](C[C@]1(C)O2)C)C")
    if not mol.HasSubstructMatch(iridoid_pattern):
        return False, "Core iridoid scaffold not found"
    
    # Check for common iridoid modifications
    
    # Glycosidic substituents
    sugar_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    
    # Acyl substituents (e.g., caffeoyl, feruloyl)
    acyl_pattern = Chem.MolFromSmarts("O=C(O)C=Cc1ccc(O)cc1")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    
    # Oxidation patterns (alcohols, ketones, carboxylic acids)
    alcohol_pattern = Chem.MolFromSmarts("[OH]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    
    ketone_pattern = Chem.MolFromSmarts("C=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    
    # Ring cleavage (secoiridoids)
    seco_pattern = Chem.MolFromSmarts("[C@H]1[C@H]2[C@@H]([C@H](C[C@@]1(C)O2)C)CC")
    seco_matches = mol.GetSubstructMatches(seco_pattern)
    
    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    if len(sugar_matches) > 0 or len(acyl_matches) > 0 or \
       len(alcohol_matches) > 0 or len(ketone_matches) > 0 or \
       len(acid_matches) > 0 or len(seco_matches) > 0 or \
       (mol_wt >= 200 and mol_wt <= 500 and n_rotatable >= 3):
        return True, "Molecule contains the core iridoid scaffold and common iridoid modifications"
    
    return False, "Molecule does not meet the criteria for an iridoid monoterpenoid"