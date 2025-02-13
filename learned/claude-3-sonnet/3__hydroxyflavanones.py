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

    # Basic flavanone core pattern - more flexible to allow substitutions
    # [#6~2] means any carbon in a ring of size 2 or more
    # This pattern matches the essential C-O-C-C(=O)-C framework
    flavanone_core = Chem.MolFromSmarts("""
        [OX2]1[CH2][CH1]([#6;R1])[CX3](=[OX1])[#6;R1]2[#6;R1]=,:[#6;R1]=,:[#6;R1]=,:[#6;R1]=,:[#6;R1]12
    """)
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "Missing flavanone core structure"

    # Pattern for 3'-hydroxyl group on the phenyl ring
    # The pattern looks for:
    # 1. The flavanone core with the 2-position carbon as the attachment point
    # 2. A phenyl ring attached at position 2
    # 3. A hydroxyl group at the meta position (3') relative to the attachment point
    hydroxyflavanone_pattern = Chem.MolFromSmarts("""
        [OX2]1[CH2][CH1]([#6;R1]2=CC(O)=CC=C2)[CX3](=[OX1])[#6;R1]3[#6;R1]=,:[#6;R1]=,:[#6;R1]=,:[#6;R1]=,:[#6;R1]13
    """)
    
    if not mol.HasSubstructMatch(hydroxyflavanone_pattern):
        # Try alternative pattern that's more flexible with substitutions
        alt_pattern = Chem.MolFromSmarts("""
            [OX2]1[CH2][CH1]([#6;R1]2[#6;R1]=,:[#6;R1](O)[#6;R1]=,:[#6;R1]=,:[#6;R1]2)[CX3](=[OX1])[#6;R1]3[#6;R1]=,:[#6;R1]=,:[#6;R1]=,:[#6;R1]=,:[#6;R1]13
        """)
        if not mol.HasSubstructMatch(alt_pattern):
            return False, "Missing 3'-hydroxy group on phenyl ring"

    # Verify it's a flavanone (not a flavone) by checking for sp3 carbon at position 2
    sp3_check = Chem.MolFromSmarts("[OX2]1[CH2][CH1][CX3](=[OX1])")
    if not mol.HasSubstructMatch(sp3_check):
        return False, "Not a flavanone (missing sp3 carbon at position 2)"

    return True, "Contains flavanone core with 3'-hydroxy substituent"