"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: polyamine
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.

    An amino group is a nitrogen atom bonded to at least one hydrogen (primary or secondary amine),
    excluding those in amides, nitro groups, nitriles, or quaternary ammonium ions.

    This function also excludes amino acids and peptides by checking for carboxylic acid groups (COOH)
    and peptide bonds (amide linkages between amino acids).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for amino groups (primary or secondary amines)
    # Exclude amides, imines, nitriles, nitro groups, and quaternary ammonium
    amino_group_smarts = "[#7;!$(N[*]=[*]);!$(N#*);!$(N*=*);!$(N(=O)(O));!$(N(~O));H1,H2]"  # N with 1 or 2 H, excluding certain groups
    amino_group_pattern = Chem.MolFromSmarts(amino_group_smarts)
    amino_groups = mol.GetSubstructMatches(amino_group_pattern)
    amino_group_count = len(amino_groups)

    # Exclude molecules with carboxylic acid groups (amino acids)
    carboxylic_acid_smarts = "[CX3](=O)[OX1H0-,OX2H1]"
    carboxylic_acid_pattern = Chem.MolFromSmarts(carboxylic_acid_smarts)
    has_carboxylic_acid = mol.HasSubstructMatch(carboxylic_acid_pattern)

    # Exclude peptides (amide bonds between amino acids)
    peptide_bond_smarts = "[NX3][CX3](=O)[NX3]"
    peptide_bond_pattern = Chem.MolFromSmarts(peptide_bond_smarts)
    has_peptide_bond = mol.HasSubstructMatch(peptide_bond_pattern)

    if has_carboxylic_acid or has_peptide_bond:
        return False, "Molecule is an amino acid or peptide"

    if amino_group_count >= 2:
        return True, f"Contains {amino_group_count} amino groups"
    else:
        return False, f"Contains {amino_group_count} amino group(s), which is less than 2"