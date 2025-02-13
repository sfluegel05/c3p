"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    An alkaloid is a naturally occurring, basic nitrogen compound (mostly heterocyclic) 
    occurring mostly in the plant kingdom, but also found in bacteria, fungi, and animals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains nitrogen
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found"

    # Check if nitrogen is part of a heterocyclic ring
    heterocyclic_nitrogen = False
    for atom in nitrogen_atoms:
        if atom.IsInRing():
            heterocyclic_nitrogen = True
            break
    if not heterocyclic_nitrogen:
        return False, "No nitrogen in heterocyclic ring"

    # Exclude molecules with exocyclic nitrogen (amines)
    for atom in nitrogen_atoms:
        if not atom.IsInRing():
            return False, "Exocyclic nitrogen (amine) detected"

    # Exclude amino acids, peptides, proteins, nucleotides, nucleic acids, amino sugars, and antibiotics
    # These are complex checks, but we can use some simple heuristics
    # For example, amino acids have both amino and carboxyl groups
    amino_acid_pattern = Chem.MolFromSmarts("[NH2]-[CX4]-[C](=O)[OH]")
    if mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Amino acid detected"

    # Peptides and proteins have multiple amide bonds
    peptide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3H]")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    if len(peptide_matches) > 1:
        return False, "Peptide or protein detected"

    # Nucleotides and nucleic acids have phosphate and sugar moieties
    nucleotide_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    if mol.HasSubstructMatch(nucleotide_pattern):
        return False, "Nucleotide or nucleic acid detected"

    # Amino sugars have both amino and hydroxyl groups on a sugar ring
    amino_sugar_pattern = Chem.MolFromSmarts("[NH2]-[CX4]-[OX2H]")
    if mol.HasSubstructMatch(amino_sugar_pattern):
        return False, "Amino sugar detected"

    # Antibiotics often have complex structures, but we can check for common antibiotic substructures
    # This is a simplified check and may not cover all cases
    antibiotic_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3H]")
    if mol.HasSubstructMatch(antibiotic_pattern):
        return False, "Antibiotic detected"

    # If all checks pass, classify as alkaloid
    return True, "Contains nitrogen in a heterocyclic ring and meets alkaloid criteria"