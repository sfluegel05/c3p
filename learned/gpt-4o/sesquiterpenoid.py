"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is any terpenoid derived from a sesquiterpene, often rearranged or modified.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a sesquiterpenoid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (10 <= c_count <= 20):
        return False, f"Contains {c_count} carbons, outside typical range for sesquiterpenoids"

    # Check for at least one ring (common in sesquiterpenoids)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, "Lacks ring structures typically present in sesquiterpenoids"

    # Check for typical sesquiterpene motifs
    sesquiterpene_motif = Chem.MolFromSmarts('C1=CC=C1')  # Replace with a more specific motif if needed
    if not mol.HasSubstructMatch(sesquiterpene_motif):
        return False, "Does not contain common sesquiterpene motif"

    return True, "Contains characteristics typical of a sesquiterpenoid"

# Example
print(is_sesquiterpenoid('O=C1C=C2C=CC(=O)[C@@]([C@]2(C)C[C@]1(O)C(=C)COC(=O)C)(O)C'))  # Example test