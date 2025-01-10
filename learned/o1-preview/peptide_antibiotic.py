"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: peptide antibiotic
"""
from rdkit import Chem

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    A peptide antibiotic is a peptide composed of amino acid residues linked via peptide bonds,
    often cyclic, and may contain unusual amino acids or modifications such as sulfur-containing rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define peptide bond pattern (amide bond between amino acids)
    peptide_bond_pattern = Chem.MolFromSmarts("N[C;!R][C;!R](=O)")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bonds)
    
    # Check for cyclic peptide (macrocycles involving peptide bonds)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    has_macrocycle = False
    for ring in atom_rings:
        if len(ring) > 8:
            # Check if ring contains peptide bonds
            ring_atoms = set(ring)
            for bond in peptide_bonds:
                if bond[0] in ring_atoms or bond[1] in ring_atoms:
                    has_macrocycle = True
                    break
            if has_macrocycle:
                break

    # Define patterns for sulfur-containing rings (e.g., thiazole, thiazoline)
    sulfur_ring_pattern = Chem.MolFromSmarts("[#16;r]")
    sulfur_rings = mol.GetSubstructMatches(sulfur_ring_pattern)
    has_sulfur_ring = len(sulfur_rings) > 0

    # Classification logic
    if num_peptide_bonds >= 6 and has_macrocycle:
        return True, f"Molecule is a cyclic peptide with {num_peptide_bonds} peptide bonds"
    elif num_peptide_bonds >= 4 and has_sulfur_ring:
        return True, f"Molecule is a peptide with sulfur-containing rings and {num_peptide_bonds} peptide bonds"
    elif num_peptide_bonds >= 8:
        return True, f"Molecule is a peptide with {num_peptide_bonds} peptide bonds"
    else:
        return False, f"Does not meet criteria: {num_peptide_bonds} peptide bonds, cyclic: {has_macrocycle}, sulfur ring: {has_sulfur_ring}"

__metadata__ = {
    'chemical_class': {
        'name': 'peptide antibiotic',
        'definition': 'A chemically diverse class of peptides that exhibit antimicrobial properties.'
    }
}