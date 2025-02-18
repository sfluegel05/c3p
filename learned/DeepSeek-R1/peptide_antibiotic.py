"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Criteria include multiple amide bonds in peptide backbone structure and molecular features.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all amide bonds (peptide bonds) using SMARTS pattern
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3H0,H1]")
    amide_count = len(mol.GetSubstructMatches(amide_pattern))
    if amide_count < 4:  # Reduced threshold to 4 for flexibility
        return False, f"Insufficient amide bonds ({amide_count} < 4)"
    
    # Check for alpha-amino acid patterns (CH-R-C(=O)-N)
    alpha_aa_pattern = Chem.MolFromSmarts("[CH1X4]([!H0])[CX3](=O)[NX3]")
    if not mol.HasSubstructMatch(alpha_aa_pattern):
        return False, "No alpha-amino acid residues detected"
    
    # Molecular weight check (500-3000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 3500:  # Slightly expanded range
        return False, f"Molecular weight {mol_wt:.1f} Da outside range (400-3500 Da)"
    
    # Check for modified amide groups (e.g., N-methyl)
    modified_amide = Chem.MolFromSmarts("[NX3H0]([#6])[CX3](=O)")  # Any N-alkyl amide
    if mol.HasSubstructMatch(modified_amide):
        return True, "Contains modified amide groups"
    
    # Check for cyclic structures containing amides
    # SMARTS for amide in a ring: amide group where either C or N is in a ring
    cyclic_amide = Chem.MolFromSmarts("[CX3](=O)[NX3H0,H1]@*")
    if mol.HasSubstructMatch(cyclic_amide):
        return True, "Contains cyclic structure with amide bonds"
    
    # Check for extended peptide chain (at least 3 consecutive amide-linked residues)
    # Pattern: N-C-C(=O)-N repeated
    peptide_chain = Chem.MolFromSmarts("[NX3][CX4H][CX3](=O)[NX3][CX4H][CX3](=O)[NX3]")
    if mol.HasSubstructMatch(peptide_chain):
        return True, "Contains extended peptide backbone"
    
    # Final check: sufficient nitrogen count and complexity
    nitrogen_count = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
    if nitrogen_count >= 6 and amide_count >= 4:
        return True, "High nitrogen/amide content characteristic of peptides"
    
    return False, "Does not meet peptide antibiotic criteria"