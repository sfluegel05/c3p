"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: CHEBI:15837 sphingomyelin
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    A sphingomyelin is a sphingolipid with a sphingoid base backbone that has:
    - A long-chain sphingoid base backbone (generally 18 carbons) with an amino group and two hydroxyl groups
    - An amide linkage between the amino group and a fatty acid
    - A phosphorylcholine group esterified to the terminal hydroxyl group

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to accurately count implicit hydrogens
    mol = Chem.AddHs(mol)

    # Define sphingoid base backbone pattern (long-chain amino diol)
    # - Long aliphatic chain with 2 or more carbons
    # - Central amino group (-NH-)
    # - Two hydroxyl groups (-OH), one at terminal position
    sphingoid_pattern = Chem.MolFromSmarts("[C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][N;R0][C;R0][O;R0]")
    sphingoid_matches = mol.GetSubstructMatches(sphingoid_pattern)
    if not sphingoid_matches:
        return False, "No sphingoid base backbone found"

    # Identify amide linkage between amino group and fatty acid
    # Pattern for amide bond attached to nitrogen
    amide_pattern = Chem.MolFromSmarts("[$([N;R0]-C(=O)-[C;R0])]")  # Non-ring nitrogen connected to carbonyl carbon
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide linkage with fatty acid found"

    # Check that amide nitrogen is part of sphingoid base
    sphingoid_atoms = set()
    for match in sphingoid_matches:
        sphingoid_atoms.update(match)
    amide_nitrogen_idx = amide_matches[0][0]
    if amide_nitrogen_idx not in sphingoid_atoms:
        return False, "Amide nitrogen not part of sphingoid backbone"

    # Identify phosphorylcholine group attached via ester linkage
    phosphorylcholine_pattern = Chem.MolFromSmarts("[O;R0]-[P](=O)([O;R0])[O;R0]-[C;R0][C;R0][N+](C)(C)C")
    phosphorylcholine_matches = mol.GetSubstructMatches(phosphorylcholine_pattern)
    if not phosphorylcholine_matches:
        return False, "No phosphorylcholine group found"

    # Check that phosphorylcholine is connected to terminal hydroxyl of sphingoid base
    phosphorylcholine_oxygen_idx = phosphorylcholine_matches[0][0]
    terminal_oxygen_idxs = [atom_idx for atom_idx in sphingoid_atoms if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 8]
    if phosphorylcholine_oxygen_idx not in terminal_oxygen_idxs:
        return False, "Phosphorylcholine group not connected to terminal hydroxyl of sphingoid base"

    return True, "Molecule is a sphingomyelin with correct structural features"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:15837',
        'name': 'sphingomyelin',
        'definition': 'Any of a class of phospholipids in which the amino group of a sphingoid base is in amide linkage with one of several fatty acids, while the terminal hydroxy group of the sphingoid base is esterified to phosphorylcholine.',
        'parents': ['CHEBI:64713', 'CHEBI:65409']
    }
}