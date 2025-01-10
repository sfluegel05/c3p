"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    Cardiac glycosides typically have a steroid skeleton with a lactone ring
    and sugar residues attached.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extended steroid backbone pattern (complex enough to capture variations)
    steroid_pattern = Chem.MolFromSmarts('C1C2CCC3C(C1)C(CC4C3C2CCCC4)C')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
     
    # Furanone lactone (common for cardiac steroids)
    lactone_pattern = Chem.MolFromSmarts('C=1OC(C=O)C=C1')  # Rings with double bonds and lactone
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Detection of glycosidic bonds (sugar attachments with O-linkages)
    sugar_pattern = Chem.MolFromSmarts('[O;R][C;R]')
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moiety found"

    # Validation with conditional logic on the oxygen count indicating complexity
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 4:
        return False, f"Found {oxygen_count} oxygens, typically too few for cardiac glycosides."
    
    return True, "Contains characteristic steroid backbone, lactone ring, and glycosidic bonds typical of cardiac glycosides."