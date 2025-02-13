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

    # Look for ester pattern (C(=O)OC)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester groups for clavulone classification, found {len(ester_matches)}"
    
    # Check for cyclopentane ring and conjugated system
    conjugated_system_pattern = Chem.MolFromSmarts("C=CC=CC")  # Simple conjugated pattern
    if mol.HasSubstructMatch(conjugated_system_pattern):
        conjugated_present = True
    else:
        conjugated_present = False
    
    # Check for halogen atoms
    halogens = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() in {9, 17, 35, 53}]
    if not halogens:
        return False, "No halogen atoms found, not typical for marine prostanoids"
    
    # Chiral centers check
    chiral_centers = Chem.FindPotentialChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 1:
        return False, f"Insufficient chiral centers, found {len(chiral_centers)}"

    reason = f"Contains ester groups, halogen atoms, and {'conjugated system' if conjugated_present else 'no conjugated system'}"
    return True, reason