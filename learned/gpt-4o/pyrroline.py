"""
Classifies: CHEBI:23763 pyrroline
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is defined as any organic heteromonocyclic compound with a structure based on a dihydropyrrole.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a five-membered ring with one nitrogen atom and one double bond
    pyrroline_pattern = Chem.MolFromSmarts("C1C=CCN1")
    pyrroline_matches = mol.GetSubstructMatches(pyrroline_pattern)
    
    if not pyrroline_matches:
        return False, "No dihydropyrrole (pyrroline) ring found"

    # Additional check for exactly one nitrogen atom in the identified ring
    for match in pyrroline_matches:
        ring = set(match)
        nitrogen_count = sum(1 for atom_idx in ring if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 7)
        if nitrogen_count == 1:
            return True, "Contains a dihydropyrrole (pyrroline) ring"

    return False, "Ring found, but does not match pyrroline specifications (e.g., wrong number of nitrogen atoms)"