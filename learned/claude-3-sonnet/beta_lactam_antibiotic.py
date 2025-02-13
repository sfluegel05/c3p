"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: CHEBI:35519 beta-lactam antibiotic
An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for beta-lactam ring pattern
    beta_lactam_pattern = Chem.MolFromSmarts("[NR]1[CR]2[CR][CR][NR]1[CR]2=O")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Look for heterocyclic rings (4-membered rings containing nitrogen)
    heterocyclic_rings = [ring for ring in mol.GetRingInfo().AtomRings()
                          if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring)]
    if not heterocyclic_rings:
        return False, "No heterocyclic rings found"

    # Look for nitrogen atoms outside the beta-lactam ring
    n_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not any(idx not in ring for ring in heterocyclic_rings for idx in ring):
        return False, "No nitrogen atoms outside the beta-lactam ring"

    return True, "Contains a beta-lactam ring and heterocyclic rings with nitrogen atoms outside the beta-lactam ring"