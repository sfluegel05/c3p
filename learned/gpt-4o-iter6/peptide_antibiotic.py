"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Peptide antibiotics are complex peptides with antimicrobial properties,
    often showing structural motifs like cyclic, branched peptides with unique
    functional groups and diverse heteroatoms.

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
    
    # Check for peptide bonds pattern: N-C(=O)
    peptide_bond_pattern = Chem.MolFromSmarts("NC(=O)")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 5:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, less than expected"
    
    # Check for common structural motifs
    extended_patterns = [
        Chem.MolFromSmarts("[C@@H]1[C@H](C)[C@H](C1)"),  # Detecting small chiral cyclic structures
        Chem.MolFromSmarts("CCC(O)=O"),  # Terminal carboxyl groups
        Chem.MolFromSmarts("C1CNC(=O)N1"),  # Beta-lactam indicative structures
        Chem.MolFromSmarts("C1CNC(=O)C1"),  # Cyclic ester indicative structures
        Chem.MolFromSmarts("C1=NC=CN=C1"),  # Pyrimidine
        Chem.MolFromSmarts("C1=CNC=CN=C1"),  # Pyridine rings
        Chem.MolFromSmarts("C1=CC=CN1"),  # Indole or pyridine rings
    ]
    
    for pattern in extended_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains known structural motifs common in peptide antibiotics"
    
    # Analyze molecular complexity and ring systems
    n_atoms = mol.GetNumAtoms()
    n_heteroatoms = sum(atom.GetAtomicNum() not in [6, 1] for atom in mol.GetAtoms())
    ring_info = mol.GetRingInfo()
    if Descriptors.MolWt(mol) > 1000 and n_heteroatoms > 10 and ring_info.NumRings() > 5:
        ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
        if any(size > 8 for size in ring_sizes) and len(set(ring_sizes)) > 1:
            return True, "Contains complex cyclic structures typical of peptide antibiotics"
    
    return False, "No clear indication of structural motifs typical of peptide antibiotics found"