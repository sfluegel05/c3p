"""
Classifies: CHEBI:28863 flavanones
"""
"""
Classifies: CHEBI:24064 flavanone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    Flavanones have a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core flavanone structure - more flexible SMARTS pattern
    # Matches the basic chroman-4-one skeleton with an aryl group at position 2
    flavanone_core = Chem.MolFromSmarts(
        "[#6]~1~2[#6](=[O:1])[#6][#6]([#6]([#6,#1])~[#6,#1;R0]~3~[#6]~[#6]~[#6]~[#6]~[#6]3)O[#6]~2=[#6]~[#6]~[#6]~1"
    )
    
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "Missing core flavanone skeleton"

    # Check for absence of double bond between positions 2 and 3
    # (this distinguishes flavanones from flavones)
    flavone_pattern = Chem.MolFromSmarts("O1c2ccccc2C(=O)C=C1c1ccccc1")
    if mol.HasSubstructMatch(flavone_pattern):
        return False, "Contains C2-C3 double bond (flavone structure)"

    # Additional validation patterns
    patterns = {
        # Check for ketone at position 4
        'ketone': Chem.MolFromSmarts("[#6]-1-[#6](=O)-[#6]-[#6]-O-[#6]-1"),
        
        # Check for aryl group at position 2 (allowing substitutions)
        'aryl': Chem.MolFromSmarts("O1[CH]([CH2]C(=O))c2[#6]~[#6]~[#6]~[#6]~[#6]2"),
        
        # Common substitution pattern (5,7-dihydroxy)
        'hydroxy': Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc2cc(O)cc(O)c21")
    }

    # Must have ketone and aryl group
    if not (mol.HasSubstructMatch(patterns['ketone']) and 
            mol.HasSubstructMatch(patterns['aryl'])):
        return False, "Missing required ketone at position 4 or aryl group at position 2"

    # Count the number of rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient ring count for flavanone structure"

    # Additional check for proper oxidation state at C2-C3
    # (should be single bonds, not double)
    c2_c3_pattern = Chem.MolFromSmarts("O1[CH]-[CH2]C(=O)")
    if not mol.HasSubstructMatch(c2_c3_pattern):
        return False, "Incorrect oxidation state at C2-C3"

    return True, "Contains 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton with correct substitution pattern"