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

    # Check for long aliphatic chain pattern
    aliphatic_chain = Chem.MolFromSmarts("C/C=C/CCCC")
    if not mol.HasSubstructMatch(aliphatic_chain):
        return False, "Missing long aliphatic chain"

    # Check for presence of ester groups
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "Missing ester group"

    # Check for presence of rings, specifically looking for aromatic/conjugated ones
    ring_info = mol.GetRingInfo()
    if not ring_info.IsAtomInRingOfSize(0, 6):
        return False, "Expected conjugated/aromatic ring not found"

    # Check for chiral centers, which are common in the examples
    chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    if chiral_centers < 1:
        return False, "Lacking required stereochemical centers"

    # Check for halogens
    halogen = Chem.MolFromSmarts("[Cl,Br,I]")
    has_halogen = mol.HasSubstructMatch(halogen)

    if has_halogen:
        return True, "Contains long aliphatic chain, ester group, cyclic structure, and halogen"
    else:
        return True, "Contains long aliphatic chain, ester group, and cyclic structure"

    return None, None