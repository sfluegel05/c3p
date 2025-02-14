"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: CHEBI:35457 monoamine

A monoamine is an aralylamino compound which contains one amino group connected to an aromatic ring by a two-carbon chain.
Monoamines are derived from aromatic amino acids like phenylalanine, tyrosine, tryptophan, and the thyroid hormones by the action of aromatic amino acid decarboxylase enzymes.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoamine(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a monoamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aromatic atoms
    aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if not aromatic_atoms:
        return False, "No aromatic ring found"

    # Look for amino group (-NH2, -NH3+, -NHR, -NR2, -NR3+)
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    if not amino_matches:
        return False, "No amino group found"

    # Check if the amino group is connected to an aromatic ring via a carbon chain
    for atom_idx in amino_matches[0]:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                continue
            path = Chem.FindAllPathsOfLengthN(mol, atom_idx, neighbor.GetIdx(), 6, useBonds=False)
            for path in path:
                if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in path[1:-1]):
                    return True, "Monoamine structure detected"

    return False, "Amino group not connected to aromatic ring via carbon chain"