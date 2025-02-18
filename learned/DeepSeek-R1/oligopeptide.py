"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a peptide containing 2-20 amino acids linked by peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find peptide bonds (amide between alpha carbons)
    peptide_bond = Chem.MolFromSmarts("[CX3](=O)-[NX3]-[CX4H]")
    peptide_matches = mol.GetSubstructMatches(peptide_bond)
    n_peptide_bonds = len(peptide_matches)
    
    # Calculate possible amino acid counts (linear: n+1, cyclic: n)
    min_aa = max(n_peptide_bonds, 2)  # Cyclic case can't have less than 2
    max_aa = n_peptide_bonds + 1
    
    if not (2 <= min_aa <= 20 or 2 <= max_aa <= 20):
        return False, f"Peptide bond count ({n_peptide_bonds}) implies {min_aa}-{max_aa} amino acids (outside 2-20 range)"

    # Check for amino acid residues (allow proline and modified termini)
    # Pattern matches [NH]-CH(R)-CO where R is any side chain
    aa_pattern = Chem.MolFromSmarts("[NX3]-[CH1;!$(C=O)]-[CX3](=O)")
    if not mol.HasSubstructMatch(aa_pattern):
        return False, "No characteristic amino acid residues found"

    # Check molecular weight (typical oligopeptide <4000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 4000:
        return False, f"Molecular weight ({mol_wt:.1f} Da) exceeds oligopeptide range"

    # Check for non-peptide amides (e.g., urea, sulfonamides)
    non_peptide_amide = Chem.MolFromSmarts("[NX3]-[CX3](=O)[!#6]")  # Amide not connected to carbon
    if mol.HasSubstructMatch(non_peptide_amide):
        return False, "Contains non-peptide amide groups"

    # If passed all checks
    aa_count = f"{min_aa}-{max_aa}" if min_aa != max_aa else str(min_aa)
    return True, f"Contains {aa_count} amino acids linked by peptide bonds"