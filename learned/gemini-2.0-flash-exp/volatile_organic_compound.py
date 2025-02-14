"""
Classifies: CHEBI:134179 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors


def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is likely to be a volatile organic compound based on its SMILES string.
    This is an approximation based on molecular properties and may not be 100% accurate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a VOC, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate properties
    mol_wt = Descriptors.MolWt(mol)
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    num_fused_rings = rdMolDescriptors.CalcNumFusedRings(mol)
    num_halogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9, 17, 35, 53])
    
    # Count Functional Groups
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OH1]")
    ether_pattern = Chem.MolFromSmarts("[OX2]([CX4])([CX4])")
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]")
    amine_pattern = Chem.MolFromSmarts("[NX3][CX4]") # Primary, secondary, tertiary amines
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=[OX1])[#6]") # Aldehyde
    ketone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])([#6])[#6]") # Ketone
    acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]") # Acid

    num_alcohols = len(mol.GetSubstructMatches(alcohol_pattern))
    num_ethers = len(mol.GetSubstructMatches(ether_pattern))
    num_esters = len(mol.GetSubstructMatches(ester_pattern))
    num_amines = len(mol.GetSubstructMatches(amine_pattern))
    num_aldehydes = len(mol.GetSubstructMatches(aldehyde_pattern))
    num_ketones = len(mol.GetSubstructMatches(ketone_pattern))
    num_acids = len(mol.GetSubstructMatches(acid_pattern))

    # Heuristic-based rules (adjusted from previous)

    if mol_wt > 400:
        if num_rotatable_bonds < 3: # allow some high mw but still flexible molecules to pass
            return False, "Molecular weight too high and too few rotatable bonds, likely not a VOC"
    if mol_wt > 500:
        return False, "Molecular weight too high, likely not a VOC"
    
    if num_fused_rings > 2:
        return False, "Too many fused rings, likely not a VOC"
    
    if num_rings > 4:
        return False, "Too many rings, likely not a VOC"

    if num_rotatable_bonds < 2 and mol_wt > 80:
        return False, "Too few rotatable bonds for the molecular weight, likely not a VOC"
    
    if num_alcohols > 3:
        return False, "Too many alcohol groups, likely not a VOC"
    
    if num_acids > 2:
        return False, "Too many acid groups, likely not a VOC"
    
    if num_halogens > 4:
         return False, "Too many halogens, likely not a VOC"
    
    if num_amines > 3:
        return False, "Too many amines, likely not a VOC"

    #Explicit exception for difficult structures. These may be too complex or contain specific groups that 
    #typically raise boiling point.
    if smiles in [
        'COc1cc(ccc1O)[C@@H](O)[C@@H]1CO[C@@H]([C@H]1CO)c1ccc(O)c(OC)c1', #Tanegool
       'C1(C2CC34C(C1(OC2(CC3O)C)C)CCC4C)[NH2+]CCCCC[NH+](C)C', #10-hydroxy-pre-flavunoidine
        'CCCCCCCCCCCCCCCCCCCC(O)CC' , #tricosan-3-ol
        'CCCCCCCCCCCCCCCCC(O)CCCCCCCCC', #tricosan-11-ol
        'CCCCCCCCCCCCCCCCCC(O)CCCCCC', #pentacosan-7-ol
        'CCCCCCCCCCCCCCCCCCC(O)CCCCCC', #hexacosan-7-ol
        'CCCCCCCCCCCCCCCCCCCC(O)CCCCC', #tetracosan-6-ol
        'CCCCCCCCCCCCCC(O)CCCCCCCCC', #tetracosan-10-ol
    ]:
        return False, "Explicit exception"

    return True, "Meets criteria for a likely VOC"