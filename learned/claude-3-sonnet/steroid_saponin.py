"""
Classifies: CHEBI:61655 steroid saponin
"""
"""
Classifies: steroid saponin
Definition: Any saponin derived from a hydroxysteroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for steroid core (four fused rings)
    steroid_pattern = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core structure found"

    # Look for glycosidic linkage pattern (C-O-C with OH groups nearby)
    glycoside_pattern = Chem.MolFromSmarts("[#6]-[OX2]-[#6;R](-[OX2H1])-[#6;R](-[OX2H1])")
    if not mol.HasSubstructMatches(glycoside_pattern):
        return False, "No glycosidic linkage found"

    # Count number of hydroxyl groups
    oh_pattern = Chem.MolFromSmarts("[OX2H1]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 2:
        return False, "Insufficient hydroxyl groups"

    # Count oxygen atoms (saponins typically have many oxygens)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Too few oxygen atoms for a saponin"

    # Check molecular weight (steroid saponins are typically large molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for steroid saponin"

    # Count rings (steroid saponins have multiple rings)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 4:
        return False, "Too few rings for steroid saponin"

    # Look for sugar moiety patterns
    pyranose_pattern = Chem.MolFromSmarts("[#6;R6]-1-[#6;R6]-[#6;R6]-[#6;R6]-[#6;R6]-[#8;R6]-1")
    if not mol.HasSubstructMatch(pyranose_pattern):
        return False, "No sugar moiety found"

    # Count carbons (steroid saponins typically have many carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for steroid saponin"

    return True, "Contains steroid core with glycosidic linkages and hydroxyl groups"