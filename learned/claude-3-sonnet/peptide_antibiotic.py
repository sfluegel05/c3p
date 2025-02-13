"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: CHEBI:35492 peptide antibiotic
A chemically diverse class of peptides that exhibit antimicrobial properties.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, MolFromSmiles
from rdkit.Chem.MolStandardize import rdMolStandardize
from typing import Tuple

def is_peptide_antibiotic(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and sanitize molecule
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Standardize molecule
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)
    
    # Check for peptide bonds (-CO-NH-)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)NCC")
    if not mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "No peptide bonds found"
    
    # Check for amino acid residues
    aa_patterns = [Chem.MolFromSmarts(p) for p in ['NC(C)C(=O)O', 'NC(CC(N)=O)C(=O)O', 'NC(CCC(=O)O)C(=O)O',
                                                   'NC(Cc1cnc[nH]1)C(=O)O', 'NC(Cc1c[nH]c2ccccc12)C(=O)O',
                                                   'NC(Cc1ccccc1)C(=O)O']]
    has_aa = any(mol.HasSubstructMatch(p) for p in aa_patterns)
    if not has_aa:
        return False, "No amino acid residues found"
    
    # Check molecular weight - peptides typically 500-5000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 5000:
        return False, f"Molecular weight {mol_wt:.0f} Da outside typical range for peptides"
    
    # Check for ring systems (common in non-ribosomal peptides)
    rings = mol.GetRingInfo().AtomRings()
    if rings:
        return True, "Contains peptide bonds, amino acid residues, and ring systems (likely non-ribosomal peptide)"
    else:
        return True, "Contains peptide bonds and amino acid residues (likely ribosomal peptide)"