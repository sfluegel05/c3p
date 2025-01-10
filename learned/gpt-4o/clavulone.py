"""
Classifies: CHEBI:36092 clavulone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    Clavulones are a class of esterified prostanoids obtained from marine corals.

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

    # Look for ester pattern (O=C-O)
    ester_pattern = Chem.MolFromSmarts("O=C-O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"
    
    # Look for cyclopentane prostanoid core pattern
    prostanoid_core_pattern = Chem.MolFromSmarts("C1CCCC1")
    if not mol.HasSubstructMatch(prostanoid_core_pattern):
        return False, "No cyclopentane prostanoid core found"
    
    # Optional: Check for halogens (indicating marine source)
    halogens = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() in {9, 17, 35, 53}]
    
    # Count the esters, cyclopentane ring and check for halogens
    if len(ester_matches) >= 2:
        reason = f"Contains cyclopentane core and {len(ester_matches)} ester groups"
        if halogens:
            reason += f" with halogen atoms present: {','.join(str(hal) for hal in halogens)}"
        return True, reason

    return False, "Insufficient ester groups for clavulone classification"