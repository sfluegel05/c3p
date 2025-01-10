"""
Classifies: CHEBI:36092 clavulone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a clavulone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for complex unsaturated chain (approximate)
    unsaturated_chain = Chem.MolFromSmarts("C=CCCC")
    if not mol.HasSubstructMatch(unsaturated_chain):
        return False, "Missing complex unsaturated carbon chain"

    # Check for cyclic structures
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings() > 0:
        return False, "Lacking necessary cyclic structures"

    # Check for ester groups with specific position
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, "Insufficient ester groups in the expected positions"

    # Check for chiral centers
    chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    if chiral_centers < 2:
        return False, "Lacking required stereochemical centers"

    # Check for presence of halogens
    halogen = Chem.MolFromSmarts("[Cl,Br,I]")
    has_halogen = mol.HasSubstructMatch(halogen)

    # Classify based on the structural features
    if has_halogen:
        return True, "Contains complex unsaturated chain, ester groups, cyclic structures, and halogen"
    else:
        return False, "Missing required halogen"

    return None, None