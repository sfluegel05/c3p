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

    # First check for the basic isoflavone core
    # This pattern represents the core chromone structure with a carbon at position 3
    isoflavone_core = Chem.MolFromSmarts("O=C1c2ccccc2OCC1")
    if not mol.HasSubstructMatch(isoflavone_core):
        return False, "No isoflavone core structure found"

    # Check for the complete isoflavone structure including the phenyl at position 3
    # The pattern includes the chromone core with a phenyl ring attached at position 3
    complete_isoflavone = Chem.MolFromSmarts("O=C1c2ccccc2OC(c3ccccc3)C1")
    if not mol.HasSubstructMatch(complete_isoflavone):
        return False, "Missing required phenyl ring at position 3"

    # Check specifically for 7-hydroxy group
    # This pattern maps the hydroxy group at position 7 of the chromone core
    hydroxy_7_pattern = Chem.MolFromSmarts("O=C1c2c(O)cccc2OC(c3ccccc3)=C1")
    if not mol.HasSubstructMatch(hydroxy_7_pattern):
        return False, "No hydroxy group at position 7"

    # Additional validation
    # Count oxygens (minimum 3: chromone O, C=O, and 7-OH)
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygens < 3:
        return False, "Insufficient number of oxygen atoms"

    # Verify the presence of two aromatic rings (A and B rings)
    num_aromatic_rings = len(mol.GetAromaticRings())
    if num_aromatic_rings < 2:
        return False, "Insufficient number of aromatic rings"

    # Check for proper conjugation in the core structure
    conjugated_pattern = Chem.MolFromSmarts("O=C1c2ccccc2OC=C1")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Missing proper conjugation in core structure"

    return True, "Valid 7-hydroxyisoflavone structure found with hydroxy group at position 7"