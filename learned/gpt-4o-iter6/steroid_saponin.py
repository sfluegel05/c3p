"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is characterized by a steroid backbone with attached sugar moieties
    and derived from a hydroxysteroid.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a steroid saponin, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Enhanced SMARTS pattern for the steroid backbone (including variations)
    steroid_pattern = Chem.MolFromSmarts("[#6]1-[#6](-[#6]2-[#6](-[#6](-[#6]3-[#6]-[#6]-[#6](-[#6]-3)-[#6]4-[#6]-[#6]-[#6](-[#6]-1)-[#6]-4)-[#6]-2)-[#6])")
    
    # Verify steroid backbone presence
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # SMARTS pattern for sugars, targeting beta-D-glucopyranoside-type connections
    sugar_pattern = Chem.MolFromSmarts("[OX2H][C@H]1[C@@H]([OX2H])[C@@H]([OX2H])[C@H]([OX2H])[C@H]1[OX2H]")
    
    # Finding attached sugar moieties
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No sugar moieties attached"
    
    # Check for hydroxyl groups on the steroid core - various positions
    hydroxy_pattern = Chem.MolFromSmarts("C(C)(C)O")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Steroid core lacks hydroxyl groups"
    
    return True, "Contains a hydroxysteroid backbone with attached sugar moieties consistent with a steroid saponin"