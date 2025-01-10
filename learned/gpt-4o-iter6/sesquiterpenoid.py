"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is derived from a sesquiterpene and typically has a C15 base
    structure possibly modified, with ring structures formed by isoprene units.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Sesquiterpenoids typically contain 15 carbons, allow slight flexibility
    if c_count < 13 or c_count > 18:
        return False, f"Expected a sesquiterpenoid-like carbon count around 15, got {c_count}"
    
    # Look for multiple isoprene units: C=C-C-C
    isoprene_unit = Chem.MolFromSmarts('C=C-C-C')
    isoprene_matches = mol.GetSubstructMatches(isoprene_unit)
    
    # Validate at least 2 isoprene matches; flexibility for rearrangements
    if len(isoprene_matches) < 2:
        return False, "Insufficient isoprene units detected"
    
    # Coarse evaluation of sesquiterpenoid substructure: general skepticism without matching core
    sesquiterpenoid_cores = [
        Chem.MolFromSmarts('C1=C[C@@H]2CC[C@H]2CC1'),  # bicyclic core with chirality
        Chem.MolFromSmarts('C1(CCCCC1)C'),             # monocyclic sesquiterpenes
        Chem.MolFromSmarts('C1[C@H]2C=C[C@H]2[C@@H]C1'), # tri- and polycyclic variations
    ]
    
    # Check for matches with major sesquiterpenoid-like cores
    core_matches = any(mol.HasSubstructMatch(core) for core in sesquiterpenoid_cores)
    if core_matches and len(isoprene_matches) >= 2:
        return True, "Contains a recognizable sesquiterpenoid core and sufficient isoprene units"
    
    return False, "Does not exhibit a recognizably sesquiterpenoid structural backbone"