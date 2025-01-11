"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    Oligopeptides generally contain a small number of amino acids, typically 2-20.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify peptide bonds (typical amide bonds in peptides)
    peptide_bond_pattern = Chem.MolFromSmarts("N-[C]-C(=O)")
    num_peptide_bonds = len(mol.GetSubstructMatches(peptide_bond_pattern))
    
    if num_peptide_bonds < 1:
        return False, "No peptide bonds identified"
    
    # Check if the number of peptide bonds make sense for an oligopeptide
    if num_peptide_bonds > 20:
        return False, f"Too many peptide bonds for an oligopeptide: {num_peptide_bonds}"
    
    # Recognize a wider variety of amino acid residues
    amino_acid_patterns = [
        Chem.MolFromSmarts("[C@@H](N)C(=O)O"),  # Typical L-alpha amino acids
        Chem.MolFromSmarts("[C@H](N)C(=O)O"),   # D-alpha amino acids
        Chem.MolFromSmarts("N[C@@]([C]=O)([C@H]C)"),  # Beta amino acids
        Chem.MolFromSmarts("NC(=O)[C@H]")      # Flattened general amide linkage
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in amino_acid_patterns):
        return False, "Missing recognizable amino acid residues"
    
    # Ensure no small inappropriate ring sizes 
    ring_info = mol.GetRingInfo()
    for atom in mol.GetAtoms():
        for ring_size in range(3, 5):  # Very small rings, rare in peptides
            if ring_info.IsAtomInRingOfSize(atom.GetIdx(), ring_size):
                return False, f"Detected a small ring of size {ring_size}, unlikely in oligopeptides"

    # Check molecular weight - oligopeptides are usually between 100-3000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 3000:
        return False, f"Molecular weight {mol_wt:.1f} Da out of range for oligopeptide"
    
    return True, "Oligopeptide structure identified with recognized peptide bonds and amino acid patterns"