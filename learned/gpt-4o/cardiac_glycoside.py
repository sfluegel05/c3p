"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    Cardiac glycosides typically have a steroid skeleton with a lactone ring
    and sugar residues.
    
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
    
    # More general steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C3C(CCC4C3CCC4C2C1)C')  # Simplified motif
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Adaptable lactone ring pattern (more generalized)
    lactone_patterns = [
        Chem.MolFromSmarts('O=C1CCCO1'),  # Simplest lactone ring
        Chem.MolFromSmarts('C1=COC(=O)C1'),  # Common furan (5-membered) lactone
        Chem.MolFromSmarts('C1=CCC(=O)O1')  # Variation of lactone structures
    ]
    lactone_found = any(mol.HasSubstructMatch(pattern) for pattern in lactone_patterns)
    if not lactone_found:
        return False, "No lactone ring found"
    
    # Detect sugar moieties (allow variability in -O- linkages)
    sugar_pattern = Chem.MolFromSmarts('[C,O][C,O].[C,O][C,O]')
    sugar_attachments = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_attachments) < 1:
        return False, "Insufficient sugar moieties found"
  
    # Further ensure that it's a cardiac glycoside by confirming number of oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 4:
        return False, f"Found {oxygen_count} oxygens, which is typically low for cardiac glycosides."
    
    return True, "Contains characteristic steroid backbone, lactone ring, and glycosidic bonds typical of cardiac glycosides."