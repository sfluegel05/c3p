"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: CHEBI:35510 macrolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolring

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is a macrocyclic lactone with a ring of twelve or more members derived from a polyketide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find largest ring size
    largest_ring_size = max([len(ring) for ring in rdmolring.GetRingInfo(mol).AtomRings()])
    if largest_ring_size < 12:
        return False, "No ring with 12 or more atoms found"

    # Check if largest ring contains a lactone
    lactone_pattern = Chem.MolFromSmarts("[O-]C(=O)")
    largest_ring_atoms = set(atom.GetIdx() for ring in rdmolring.GetRingInfo(mol).AtomRings() if len(ring) == largest_ring_size for atom in ring)
    if not any(mol.HasSubstructMatch(lactone_pattern, atoms=list(largest_ring_atoms))):
        return False, "Largest ring does not contain a lactone"

    # Check if molecule is derived from a polyketide
    polyketide_pattern = Chem.MolFromSmarts("[C](-[C])(-[C])(-[C])-[C]")
    if not mol.HasSubstructMatch(polyketide_pattern):
        return False, "Not derived from a polyketide"

    return True, "Macrocyclic lactone with a ring of twelve or more members derived from a polyketide"