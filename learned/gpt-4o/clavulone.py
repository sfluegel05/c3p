"""
Classifies: CHEBI:36092 clavulone
"""
from rdkit import Chem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    Clavulones are esterified prostanoids from marine corals, often containing
    ester groups, halogens, and conjugated systems.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester pattern (C(=O)OC) - Clavulones might require more than one ester group
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester groups for clavulone classification, found {len(ester_matches)}"

    # Check for broader conjugated systems
    conjugated_system_pattern = Chem.MolFromSmarts("C=C-,-{2,}")
    if not mol.HasSubstructMatch(conjugated_system_pattern):
        return False, "Lacks complex conjugated system typical of clavulones"
    
    # Check for halogen atoms (Cl, Br, I)
    halogens = any(atom.GetAtomicNum() in {9, 17, 35, 53} for atom in mol.GetAtoms())
    if not halogens:
        return False, "No halogen atoms found, not typical for marine prostanoids"
    
    # Check for five-membered rings
    ring_info = mol.GetRingInfo()
    five_membered_rings = any(len(ring) == 5 for ring in ring_info.AtomRings())
    if not five_membered_rings:
        return False, "No five-membered rings found, atypical structure for prostanoids"
    
    # Check for chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return False, "No chiral centers indicated in SMILES"

    # If all checks pass, classify as clavulone
    reason = "Contains ester groups, halogen atoms, complex conjugated system, and five-membered rings with chiral centers"
    return True, reason