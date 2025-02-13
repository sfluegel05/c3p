"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: CHEBI:16670 oligopeptide
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

AA_CODES = "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a peptide containing a relatively small number of amino acids.

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

    # Check molecular weight range (300-2000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 2000:
        return False, f"Molecular weight ({mol_wt:.2f} Da) outside typical oligopeptide range"

    # Count amino acid residues
    aa_residues = sum(1 for atom in mol.GetAtoms() if atom.GetSmarts()[0] in AA_CODES)
    if aa_residues < 2:
        return False, "Less than 2 amino acid residues"

    # Look for peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("[N;X3;H2,H1][C;X4][C;X3](=O)[N;X3]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 1:
        return False, "No peptide bonds found"

    # Check for common amino acid sidechains
    sidechain_patterns = [Chem.MolFromSmarts(smi) for smi in [
        "[NH2]-[CH2]-[CH2]-[CH2]-[CH2]-[NH]",  # Lys
        "[CH3]-[CH2]-[CH]-[NH]",  # Ala
        "[CH3]-[CH]-[CH3]",  # Val
        "[CH3]-[CH2]-[CH2]-[CH]-[NH]",  # Leu
        "[CH2]-[CH2]-[CH2]-[CH2]-[NH]",  # Pro (cyclic)
        "[CH3]-[S]-[CH3]",  # Met
        "[cH]-1-[cH]-[cH]-[cH]-[cH]-[cH]-1-[NH]",  # Phe, Tyr, Trp
        "[NH]-[CH2]-[CH2]-[CH2]-[NH2]",  # Arg
        "[CH2]-[CH2]-[C]-[NH]",  # Asp, Glu
        "[CH2]-[C]-[NH2]",  # Asn, Gln
        "[CH2]-[CH2]-[CH2]-[NH2]",  # His
        "[S]-[CH3]"  # Cys
    ]]
    has_sidechains = any(mol.HasSubstructMatch(pattern) for pattern in sidechain_patterns)
    if not has_sidechains:
        return False, "No common amino acid sidechains found"

    return True, "Contains multiple amino acid residues connected by peptide bonds"