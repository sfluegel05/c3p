"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: 7-hydroxyisoflavones
A hydroxyisoflavone compound having a hydroxy group at the 7-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic isoflavone core pattern:
    # - Benzopyrone system with a carbon attached at position 3
    # - More flexible pattern that allows for substitutions
    isoflavone_core = Chem.MolFromSmarts("[#6]1=[#6]-[#6](=O)-c2c(o1)cccc2")
    if not mol.HasSubstructMatch(isoflavone_core):
        return False, "No isoflavone core structure found"

    # Check for aromatic ring at position 3
    # More flexible pattern that allows for substituted phenyl rings
    c3_aromatic = Chem.MolFromSmarts("[#6]1=[#6]-[#6](=O)-c2c(o1)cccc2-[#6]3:[#6]:[#6]:[#6]:[#6]:[#6]:3")
    if not mol.HasSubstructMatch(c3_aromatic):
        return False, "Missing required aromatic ring at position 3"

    # Pattern for 7-hydroxy position
    # More flexible pattern that accounts for the core structure with OH at position 7
    # The pattern maps the whole core structure to ensure correct position
    hydroxy_7_pattern = Chem.MolFromSmarts("O-[c]1[c][c][c]2[c]([c]1)[o][c][c](-[#6]:3:[#6]:[#6]:[#6]:[#6]:[#6]:3)[c]2=O")
    if not mol.HasSubstructMatch(hydroxy_7_pattern):
        return False, "No hydroxy group at position 7"

    # Additional validation checks
    # Count oxygens (minimum 3: chromone O, C=O, and 7-OH)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygens < 3:
        return False, "Insufficient number of oxygen atoms"

    # Verify aromatic nature (minimum 10 aromatic atoms for basic isoflavone)
    num_aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if num_aromatic_atoms < 10:
        return False, "Insufficient aromatic system"

    return True, "Valid 7-hydroxyisoflavone structure found"