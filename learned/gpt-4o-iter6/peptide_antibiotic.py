"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Peptide antibiotics are complex peptides with antimicrobial properties,
    often showing structural motifs like cyclic or branched peptides with
    unique functional groups and diverse heteroatoms.

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
    
    # Check for peptide bond pattern: N-C(=O)
    peptide_bond_pattern = Chem.MolFromSmarts("NC(=O)")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 5:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, less than expected for a complex peptide antibiotic"
    
    # Estimates molecular complexity via atom types and counts
    n_atoms = mol.GetNumAtoms()
    n_heteroatoms = sum(atom.GetAtomicNum() != 6 for atom in mol.GetAtoms())
    
    # Typical functional groups and motifs
    unique_group_patterns = [
        Chem.MolFromSmarts("C1SCCN1"),  # Thiazoline ring
        Chem.MolFromSmarts("n1c2ccc[nH]c2c[nH]c1"),  # Indole-like
        Chem.MolFromSmarts("C1=NC=CN=C1"),  # Pyrimidine ring
        Chem.MolFromSmarts("C1=CN=CO1"),  # Oxazole/thiazole ring
    ]
    
    for pattern in unique_group_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains unique group common in peptide antibiotics"
       
    # Check for complex cyclic structure and diversity
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
        # Consider complexity with diverse ring sizes
        if len(set(ring_sizes)) > 1 and n_atoms > 30 and n_heteroatoms > 5:
            return True, "Contains complex cyclic structure and diverse elements typical of peptide antibiotics"

    return False, "No clear indication of structural motifs typical of peptide antibiotics found"