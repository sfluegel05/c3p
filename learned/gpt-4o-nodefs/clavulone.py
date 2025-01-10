"""
Classifies: CHEBI:36092 clavulone
"""
from rdkit import Chem

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

    # Common clavulone substructure: complex polyunsaturated lactone with ester groups
    clavulone_pattern = Chem.MolFromSmarts("C1(C(C(=O)C=C1)=CC)=C")
    if not mol.HasSubstructMatch(clavulone_pattern):
        return False, f"Missing characteristic polyunsaturated lactone structure"

    # Adjust code to prioritize structures found in the clavulones provided
    # Verify presence of ester groups connected to cyclic moieties and unsaturated backbones
    ester_cyclic_pattern = Chem.MolFromSmarts("C(=O)O[C@H]1C[CH]C(=O)C(=C1)C")
    if not mol.HasSubstructMatch(ester_cyclic_pattern):
        return False, f"No cyclic connected esters matching clavulone structure"

    # Check for at least one chiral center
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 1:
        return False, "Lacking required chiral center"

    # Classify based on the structural features
    return True, "Matches clavulone structure with polyunsaturated lactone core and ester groups"