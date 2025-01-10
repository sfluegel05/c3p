"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: CHEBI:35610 azole
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    An azole is a monocyclic heteroarene consisting of a five-membered ring containing nitrogen.
    The ring can also contain one or more other non-carbon atoms, such as nitrogen, sulfur, or oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the azole ring pattern: a five-membered ring with at least one nitrogen
    azole_pattern = Chem.MolFromSmarts("[n]1[c,s,o][c,s,o][c,s,o][c,s,o]1")
    if not mol.HasSubstructMatch(azole_pattern):
        return False, "No five-membered ring with at least one nitrogen found"

    # Check if the ring is monocyclic
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) != 1:
        return False, "Molecule is not monocyclic"

    # Check if the ring is aromatic
    if not mol.GetAtomWithIdx(0).GetIsAromatic():
        return False, "Ring is not aromatic"

    # Count the number of nitrogen atoms in the ring
    ring_atoms = ring_info.AtomRings()[0]
    nitrogen_count = sum(1 for atom_idx in ring_atoms if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 7)
    if nitrogen_count < 1:
        return False, "Ring does not contain any nitrogen atoms"

    # Check for other heteroatoms (sulfur or oxygen) in the ring
    heteroatom_count = sum(1 for atom_idx in ring_atoms if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() in [8, 16])
    if heteroatom_count > 0:
        return True, "Contains a five-membered aromatic ring with nitrogen and other heteroatoms"
    else:
        return True, "Contains a five-membered aromatic ring with nitrogen"

# Example usage:
# print(is_azole("O1C(=NC(=C1CC)C)CCCCCC"))  # Should return True for 5-Ethyl-2-hexyl-4-methyloxazole