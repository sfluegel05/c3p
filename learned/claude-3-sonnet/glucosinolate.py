"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: glucosinolate compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the core N-O-SO3 group (both charged and uncharged forms)
    sulfate_pattern = Chem.MolFromSmarts('[NX2]-[OX2]-[SX4](=[OX1])(=[OX1])[OX1,OX2-]')
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "Missing N-O-SO3 group"

    # Look for C=N-O linkage
    cn_bond_pattern = Chem.MolFromSmarts('[CX3]=[NX2]-[OX2]')
    if not mol.HasSubstructMatch(cn_bond_pattern):
        return False, "Missing C=N-O linkage"

    # Check for thioglycosidic linkage (C-S-C)
    thio_pattern = Chem.MolFromSmarts('[CX4]-[SX2]-[CX3]=[NX2]')
    if not mol.HasSubstructMatch(thio_pattern):
        return False, "Missing thioglycosidic linkage"

    # Check for pyranose ring with correct substitution pattern
    # More flexible pattern for beta-D-glucopyranose
    sugar_pattern = Chem.MolFromSmarts('[OX2][CX4H1][CX4H1]1[OX2][CX4H1]([SX2])[CX4H1]([OX2])[CX4H1]([OX2])[CX4H1]([OX2])[CX4H1]1[CX4H2][OX2]')
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "Missing or incorrect glucose moiety"

    # Verify the complete core structure connectivity
    core_structure = Chem.MolFromSmarts('[CX4]-[OX2]-[CX4H1]-1-[OX2]-[CX4H1](-[SX2]-[CX3](=[NX2]-[OX2]-[SX4](=[OX1])(=[OX1])[OX1,OX2-])-*)-[CX4H1](-[OX2])-[CX4H1](-[OX2])-[CX4H1](-[OX2])-[CX4H1]-1-[CX4H2]-[OX2]')
    if not mol.HasSubstructMatch(core_structure):
        return False, "Incorrect overall connectivity"

    # Count key atoms to ensure reasonable composition
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if s_count < 2:  # Need at least 2 sulfur atoms
        return False, "Insufficient sulfur atoms"
    if o_count < 7:  # Need multiple oxygens for sugar + sulfate
        return False, "Insufficient oxygen atoms"
    if n_count != 1:  # Need exactly 1 nitrogen
        return False, "Incorrect number of nitrogen atoms"

    return True, "Contains complete glucosinolate structure with correct connectivity"