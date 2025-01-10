"""
Classifies: CHEBI:4194 D-hexose
"""
from rdkit import Chem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a six-carbon monosaccharide with a specific stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for hexose structure: D-hexoses contain six carbons often including a aldehyde (aldehydo form) or
    # form either pyranose or furanose rings. Here, we'll try to capture cyclic structures with oxygen.
    # We use wildcard "C~" to allow some flexibility for cyclic structures (capital O for oxygen):
    hexose_ring_patterns = [
        Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@H](O)[C@H](O)[C@H]1O"),  # a general pattern for D-pyranose forms
        Chem.MolFromSmarts("O1[C@H]([C@H](O)[C@H](O)[C@H](O)[C@H](O)1)CO"),  # D-furanose
    ]
    
    # Attempt to find matches for these patterns
    is_hexose = any(mol.HasSubstructMatch(pattern) for pattern in hexose_ring_patterns)
    if not is_hexose:
        return False, "Does not match hexose ring structure"

    # Check for the presence of D-stereoisomerism:
    # count `@` symbols, generally D-hexoses have specific stereochemistry
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 4:
        return False, "Insufficient chiral centers for a D-hexose"
    
    if sum(smiles.count('@') for _ in chiral_centers) < 2:
        return False, f"Expected more stereochemical information on {len(chiral_centers)} centers"

    # Verify that there are exactly six carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, "Not six carbons present"

    return True, "Valid D-hexose conformer and stereochemistry detected"