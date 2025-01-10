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
    A peptide antibiotic is a peptide composed of multiple amino acid residues linked via peptide bonds.

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
    peptide_bond_pattern = Chem.MolFromSmarts("NC(=O)")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bonds)
    
    if num_peptide_bonds < 4:
        return False, f"Only {num_peptide_bonds} peptide bonds found, need at least 4"
    
    # Optionally, check for cyclic peptide (some peptide antibiotics are cyclic)
    ring_info = mol.GetRingInfo()
    has_macrocycle = any(len(ring) > 8 for ring in ring_info.AtomRings())
    if has_macrocycle:
        cyclic = "cyclic peptide"
    else:
        cyclic = "linear peptide"
    
    return True, f"Molecule is a {cyclic} with {num_peptide_bonds} peptide bonds"
    
__metadata__ = {
    'chemical_class': {
        'name': 'peptide antibiotic',
        'definition': 'A chemically diverse class of peptides that exhibit antimicrobial properties.'
    }
}