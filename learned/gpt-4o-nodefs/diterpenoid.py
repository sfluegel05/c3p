"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    Diterpenoids typically contain around 20 carbons and exhibit complex ring structures 
    and a variety of functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 18 or num_carbons > 28:
        return False, f"Contains {num_carbons} carbon atoms, most diterpenoids have carbon count close to 20, but larger range is considered"

    # Get the ring information about the molecule
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 2:
        return False, "Diterpenoids typically have multiple ring structures"

    # Check for presence of diverse functional groups
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    ketone_pattern = Chem.MolFromSmarts("C=O")
    ether_pattern = Chem.MolFromSmarts("[OX2]([#6])[#6]")
    isocyano_pattern = Chem.MolFromSmarts("[NX1]#[CX2]")
    ester_pattern = Chem.MolFromSmarts("COC(=O)")
    olefin_pattern = Chem.MolFromSmarts("C=C")  # Double bonded carbons

    has_functional_group = any([
        mol.HasSubstructMatch(alcohol_pattern),
        mol.HasSubstructMatch(ketone_pattern),
        mol.HasSubstructMatch(ether_pattern),
        mol.HasSubstructMatch(isocyano_pattern),
        mol.HasSubstructMatch(ester_pattern),
        mol.HasSubstructMatch(olefin_pattern),
    ])

    if not has_functional_group:
        return False, "Lacks common functional groups like alcohols, ketones, ethers, isocyano groups, esters, or alkenes"

    return True, "Contains characteristics typical of diterpenoids, such as complex ring structures and varied functional groups"