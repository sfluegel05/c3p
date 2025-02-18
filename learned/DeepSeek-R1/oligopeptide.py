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
    An oligopeptide is a peptide containing a relatively small number of amino acids (typically 2-20).

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

    # Check for peptide bonds (amide groups)
    peptide_bond = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    matches = mol.GetSubstructMatches(peptide_bond)
    if not matches:
        return False, "No peptide bonds found"

    # Count peptide bonds to estimate amino acid count (n_bonds = n_amino_acids - 1)
    n_amino_acids = len(matches) + 1
    if not (2 <= n_amino_acids <= 20):
        return False, f"Number of amino acids ({n_amino_acids}) outside oligopeptide range (2-20)"

    # Check for amino acid residues (at least one amino and carboxyl group per residue)
    # Using simplified check for alpha-amino acid pattern
    amino_acid_pattern = Chem.MolFromSmarts("[NH2,NH1][CH][C](=O)")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Missing characteristic amino acid structure"

    # Check molecular weight consistency (optional but adds robustness)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 2000:  # Arbitrary upper limit for oligopeptides
        return False, f"Molecular weight ({mol_wt:.1f} Da) too high for oligopeptide"

    return True, f"Contains {n_amino_acids} amino acids linked by peptide bonds"