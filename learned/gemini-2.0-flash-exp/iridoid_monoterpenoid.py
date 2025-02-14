"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): (True, reason) if molecule is an iridoid monoterpenoid,
                         (False, reason) otherwise.
                         (None, None) if error.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # More flexible Iridoid core structure (cyclopentane fused with pyran)
    # Allowing for various substitutions (X) and double bonds (=).  
    iridoid_core_pattern = Chem.MolFromSmarts("[CX4]1[CX4]([CX4])[CX4]~[OX2]~[CX4]2[CX4]1[CX4]~[CX4]2") # Basic fusion
    iridoid_core_pattern2 = Chem.MolFromSmarts("[CX4]1[CX4](=[CX4])[CX4]~[OX2]~[CX4]2[CX4]1[CX4]~[CX4]2") # with double bond

    # Secoiridoid (cleaved cyclopentane ring, pyran) - some flexibility in cleavage location
    secoiridoid_core_pattern = Chem.MolFromSmarts("[CX4]1[CX4][CX4]~[OX2]~[CX4]2[CX4]1[CX4]~[CX4]2") #sec, broken bond
    secoiridoid_core_pattern2 = Chem.MolFromSmarts("[CX4]1[CX4](=[CX4])[CX4]~[OX2]~[CX4]2[CX4]1[CX4]~[CX4]2") #sec, broken bond, one double bond

    if not (mol.HasSubstructMatch(iridoid_core_pattern) or
            mol.HasSubstructMatch(iridoid_core_pattern2) or
            mol.HasSubstructMatch(secoiridoid_core_pattern) or
            mol.HasSubstructMatch(secoiridoid_core_pattern2)):
        return False, "No core iridoid or secoiridoid structure found"
    
    # check if cyclopentane and pyran are fused or linked by a chain
    if mol.HasSubstructMatch(iridoid_core_pattern) or mol.HasSubstructMatch(iridoid_core_pattern2):
        match = mol.GetSubstructMatches(iridoid_core_pattern)
        if not match:
             match = mol.GetSubstructMatches(iridoid_core_pattern2)
        if match:
            for match_tuple in match:
                atom_indices = list(match_tuple)
                #Check if they are bonded directly
                if not mol.GetBondBetweenAtoms(atom_indices[0], atom_indices[5]): #index 0 and 5 define the fusion
                   if not mol.GetBondBetweenAtoms(atom_indices[1], atom_indices[5]): # if not a regular fusion, then a ring formed in place of a bond.
                       return False, "Core cyclopentane and pyran not fused in iridoid"
            
    if mol.HasSubstructMatch(secoiridoid_core_pattern) or mol.HasSubstructMatch(secoiridoid_core_pattern2):
        match = mol.GetSubstructMatches(secoiridoid_core_pattern)
        if not match:
            match = mol.GetSubstructMatches(secoiridoid_core_pattern2)
        if match:
            for match_tuple in match:
                atom_indices = list(match_tuple)
                 #Check if the 2 rings are linked by 1 bond
                if mol.GetBondBetweenAtoms(atom_indices[0], atom_indices[5]):
                   return False, "Core cyclopentane and pyran must be open in a seco-iridoid"
    

    return True, "Iridoid or secoiridoid core structure found"