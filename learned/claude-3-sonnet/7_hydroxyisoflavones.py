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

    # Basic chromone core pattern (benzopyrone system)
    # The c1ccc2c(c1) ensures we have the first benzene ring
    # occ(*)c2=O ensures we have the pyrone ring with something (*) at position 3
    chromone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)occ(*)c2=O")
    if not mol.HasSubstructMatch(chromone_pattern):
        return False, "No chromone core structure found"

    # Check for phenyl ring at position 3 (isoflavone specific)
    # This pattern ensures the phenyl ring is connected to the chromone at the correct position
    isoflavone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)occ(-c3ccccc3)c2=O")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Not an isoflavone - missing or incorrect phenyl ring at position 3"

    # Pattern for 7-hydroxy position
    # This pattern specifically looks for OH at position 7 of the chromone system
    hydroxy_7_pattern = Chem.MolFromSmarts("Oc1ccc2c(c1)occ(-[c])c2=O")
    if not mol.HasSubstructMatch(hydroxy_7_pattern):
        return False, "No hydroxy group at position 7"

    # Additional validation checks
    
    # Count oxygens (minimum 3: chromone O, C=O, and 7-OH)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygens < 3:
        return False, "Insufficient number of oxygen atoms"

    # Verify aromatic nature (minimum 12 aromatic atoms for basic isoflavone)
    num_aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if num_aromatic_atoms < 12:
        return False, "Insufficient aromatic system"

    # Check for required carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("O=C1c2c(occ1)cccc2")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Missing required carbonyl group"

    return True, "Valid 7-hydroxyisoflavone structure found"