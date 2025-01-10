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

    # Sesquiterpenoids typically contain around 15 carbons, allow broader flexibility
    if c_count < 12 or c_count > 20:
        return False, f"Expected a sesquiterpenoid-like carbon count around 15, got {c_count}"
    
    # Broaden definition of isoprene-related units: C=C-C or similar variations
    isoprene_variations = [
        Chem.MolFromSmarts('C=C-C'),           # basic hinge of isoprene unit
        Chem.MolFromSmarts('C=C-C-C'),         # longer units
        Chem.MolFromSmarts('C-C=C')            # inverse linkage
    ]
    
    # Validate if any isoprene-like matches occur
    isoprene_matches = any(len(mol.GetSubstructMatches(pattern)) > 0 for pattern in isoprene_variations)
    
    if not isoprene_matches:
        return False, "Insufficient isoprene-like units detected"
    
    # Expand the collection of sesquiterpenoid-like cores
    sesquiterpenoid_cores = [
        Chem.MolFromSmarts('C1=C[C@@H]2CC[C@H]2CC1'),  # bicyclic core with chirality
        Chem.MolFromSmarts('C1(CCCCC1)C'),             # monocyclic sesquiterpenes
        Chem.MolFromSmarts('C1[C@H]2C=C[C@H]2[C@@H]C1'), # tri- and polycyclic variations
        Chem.MolFromSmarts('C1=C(C)C2=C(C1)C(C2)C'),   # Additional complex cores
        Chem.MolFromSmarts('C1(CCC[C@@H]1)C')           # Cyclized simple units
    ]
    
    # Check for matches with major sesquiterpenoid-like cores
    core_matches = any(mol.HasSubstructMatch(core) for core in sesquiterpenoid_cores)
    if core_matches and isoprene_matches:
        return True, "Contains a recognizable sesquiterpenoid core and sufficient isoprene-related patterns"
    
    return False, "Does not exhibit a recognizably sesquiterpenoid structural backbone"