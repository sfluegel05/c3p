"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone is expected to possess a flavanone structure with a 4'-hydroxy group on the phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identification of flavanone core structure with attached aromatic rings
    # More generalized SMARTS for flavanone with aromatic rings (allowing substitution)
    flavanone_core_pattern = Chem.MolFromSmarts("[C@@H]1(C=O)O[C@@H]2ccc(cc2)[C@H]1")  # Generalized flavanone core
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "No flavanone core structure found"

    # Look for a 4'-hydroxy group on an aromatic ring (para position)
    hydroxy_pattern = Chem.MolFromSmarts("cc(O)cc")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "4'-hydroxy group not identified on any aromatic ring"

    # Check for specific chiral centers reflecting correct stereochemistry
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=False)
    expected_chirality = len([center for center in chiral_centers if center[1] == "1"])
    if expected_chirality < 1:
        return False, "Stereochemistry not adequately identified"

    return True, "Contains flavanone core structure with 4'-hydroxy group and correct stereochemistry"