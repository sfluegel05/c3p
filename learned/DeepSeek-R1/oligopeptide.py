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
    An oligopeptide contains 2-20 amino acids linked by peptide bonds.

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
    # Updated pattern: [CX3H1,CX4H2](=O)-[NX3H1] connected to alpha carbon (adjacent to CO and N)
    peptide_bond = Chem.MolFromSmarts("[CX3H1,CX4H2](=O)-[NX3H1]-[CX4H][!$(C=O)]")
    peptide_matches = mol.GetSubstructMatches(peptide_bond)
    n_peptide_bonds = len(peptide_matches)
    
    # Calculate possible amino acid count (n_peptide_bonds +1 for linear, n_peptide_bonds for cyclic)
    if n_peptide_bonds == 0:
        return False, "No peptide bonds found"
    
    min_aa = n_peptide_bonds  # Cyclic case
    max_aa = n_peptide_bonds + 1  # Linear case
    
    if not (2 <= min_aa <= 20 or 2 <= max_aa <= 20):
        return False, f"Peptide bonds suggest {min_aa}-{max_aa} amino acids (outside 2-20 range)"

    # Check for at least two amino acid residues with proper backbone
    # Pattern matches amino acid backbone: [NH]-CH(R)-CO (allowing proline)
    aa_pattern = Chem.MolFromSmarts("[NX3H1,NH0]-[CH1,CH2;!$(C=O)]-[CX3](=O)")
    aa_matches = len(mol.GetSubstructMatches(aa_pattern))
    if aa_matches < 2:
        return False, "Insufficient amino acid residues detected"

    # Check for non-peptide amides in main chain (exclude side chains)
    # Non-peptide amide is any amide not part of the peptide backbone pattern
    all_amides = mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)-[NX3]"))
    if len(all_amides) > n_peptide_bonds:
        return False, "Contains non-peptide amide groups"

    # Check molecular weight (typical oligopeptide <4000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 4000:
        return False, f"Molecular weight ({mol_wt:.1f} Da) exceeds oligopeptide range"

    # If passed all checks
    aa_count = f"{min_aa}-{max_aa}" if min_aa != max_aa else str(min_aa)
    return True, f"Contains {aa_count} amino acids linked by peptide bonds"