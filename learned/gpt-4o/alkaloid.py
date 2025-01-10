"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    An alkaloid is a naturally occurring, basic nitrogen compound (mostly heterocyclic).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of nitrogen atoms
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found"
    
    # Check if nitrogen is part of a heterocyclic ring
    nitrogen_in_ring = False
    for atom in nitrogen_atoms:
        if atom.IsInRing():
            nitrogen_in_ring = True
            break
    if not nitrogen_in_ring:
        return False, "Nitrogen atom is not part of a heterocyclic ring"

    # Check for exclusion cases (like amino acids, peptides)
    # Simplistic check: do not contain carboxyl group
    carboxyl_group = Chem.MolFromSmarts('C(=O)[O]')
    if mol.HasSubstructMatch(carboxyl_group):
        return False, "Contains carboxyl group, possibly an amino acid or peptide"

    return True, "Contains basic nitrogen in a heterocyclic ring, classifying as an alkaloid"