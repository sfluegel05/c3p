"""
Classifies: CHEBI:13248 anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is characterized by an aromatic ring (like a phenyl group) 
    directly connected to the nitrogen atom of an amide group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved SMARTS pattern for anilides: Ensures that an aromatic ring is connected to an amide nitrogen
    # This considers the aromatic nature of the phenyl group without specifying hydrogens
    # [c]: Aromatic any carbon, [N]: Any nitrogen, [C](=O): Carbonyl attached to nitrogen
    anilide_pattern = Chem.MolFromSmarts("c-N-C(=O)")

    match = mol.HasSubstructMatch(anilide_pattern)

    if not match:
        return False, "No anilide pattern found (aromatic ring bound to amide nitrogen)"

    return True, "Contains an aromatic ring directly bonded to amide nitrogen, fitting anilide definition"