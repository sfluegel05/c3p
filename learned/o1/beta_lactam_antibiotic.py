"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: CHEBI:35627 beta-lactam antibiotic
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic is an organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.

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

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    beta_lactam_found = False

    # Loop over atom rings to find beta-lactam ring
    for ring in atom_rings:
        if len(ring) == 4:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            num_nitrogen = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 7)
            num_carbonyl = 0
            # Check for carbonyl group (C=O) in ring
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # Carbon atom
                    for bond in atom.GetBonds():
                        if bond.GetBondType() == rdchem.BondType.DOUBLE:
                            neighbor = bond.GetOtherAtom(atom)
                            if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() in ring:
                                num_carbonyl += 1
            if num_nitrogen == 1 and num_carbonyl == 1:
                beta_lactam_found = True
                break

    if not beta_lactam_found:
        return False, "No beta-lactam ring found"

    # Check if molecule is an organonitrogen heterocycle (contains nitrogen in a ring)
    organonitrogen_heterocycle = False
    for ring in atom_rings:
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        if any(atom.GetAtomicNum() == 7 for atom in atoms_in_ring):
            organonitrogen_heterocycle = True
            break

    if not organonitrogen_heterocycle:
        return False, "Molecule is not an organonitrogen heterocycle"

    # Optional: Check for antibiotic-like properties (presence of carboxylic acid or ester group)
    antibiotic_features = False
    carboxylic_acid_smarts = '[CX3](=O)[OX1H0-,OX2H1]'
    carboxylic_acid_pattern = Chem.MolFromSmarts(carboxylic_acid_smarts)
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        antibiotic_features = True

    if not antibiotic_features:
        return False, "No carboxylic acid or ester group found"

    return True, "Contains beta-lactam ring and is an organonitrogen heterocycle"