"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: 3'-hydroxyflavanones
Compounds with a flavanone core structure and a hydroxy substituent at position 3' of the phenyl ring
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic flavanone core (chroman-4-one)
    # More flexible pattern that allows for substitutions
    flavanone_core = Chem.MolFromSmarts("O1CC(=O)c2ccccc2O1")
    if flavanone_core is None:
        return False, "Error in flavanone SMARTS pattern"
    
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "Missing flavanone core structure"

    # Check for 2-phenyl-chroman-4-one with 3'-OH
    # This pattern specifically looks for:
    # - The flavanone core
    # - A phenyl ring at position 2
    # - A hydroxyl group at the 3' position
    pattern = Chem.MolFromSmarts("""
        [O;D2]1[CH2][CH1]([#6;R1]2)C(=O)c3c(O1)cccc3
        [$([#6]:c:c:c(O):c)]
    """)
    
    if pattern is None:
        return False, "Error in SMARTS pattern"

    # Alternative pattern focusing on the 3'-OH phenyl ring
    alt_pattern = Chem.MolFromSmarts("""
        [OH]-[c]1[cH,$(c[#6])]c([cH,$(c[#6])])c([CH1]2OCC(=O)c3ccccc3O2)cc1
    """)
    
    if alt_pattern is None:
        return False, "Error in alternative SMARTS pattern"

    if mol.HasSubstructMatch(pattern) or mol.HasSubstructMatch(alt_pattern):
        # Additional check for ketone group to ensure we're not matching flavones
        ketone_pattern = Chem.MolFromSmarts("O=C1CC(Oc2ccccc12)")
        if ketone_pattern is None:
            return False, "Error in ketone SMARTS pattern"
            
        if mol.HasSubstructMatch(ketone_pattern):
            return True, "Contains flavanone core with 3'-hydroxy substituent"
    
    return False, "Missing required 3'-hydroxy flavanone structure"