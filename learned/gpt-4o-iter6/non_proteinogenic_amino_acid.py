"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determine if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    
    A non-proteinogenic amino acid is defined as one with both amino and carboxyl groups,
    plus non-standard modifications in side chains or additional groups not seen in
    the 20 standard amino acids.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a non-proteinogenic amino acid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for primary, secondary, and tertiary amines
    amino_group_pattern = Chem.MolFromSmarts('[NX3;H2,H1,H0]')
    # SMARTS pattern for carboxyl group
    carboxyl_group_pattern = Chem.MolFromSmarts('C(=O)[O-,O]')
    
    # Check for amino and carboxyl groups presence
    has_amino_group = mol.HasSubstructMatch(amino_group_pattern)
    has_carboxyl_group = mol.HasSubstructMatch(carboxyl_group_pattern)
    
    if not (has_amino_group and has_carboxyl_group):
        return False, "Must contain both amino and carboxyl groups"

    # SMARTS for oxidative modifications and added complexity in side chains
    unique_modifications = [
        Chem.MolFromSmarts('C=C'),    # Double bonds in unusual places
        Chem.MolFromSmarts('N=C'),    # Amide or imine modifications
        Chem.MolFromSmarts('N-[NX3]'), # Diverse nitrogen reduction functionality
        Chem.MolFromSmarts('S'),      # Sulfur atoms outside standard placement
        Chem.MolFromSmarts('[F,Cl,Br,I]'), # Halogenation patterns
        Chem.MolFromSmarts('O=C-N'),  # Acetyl groups or modified amidation
        Chem.MolFromSmarts('[n,N]'),  # Aromatic nitrogens
        Chem.MolFromSmarts('[CH2]O[CH2]'), # Ether linkages
        Chem.MolFromSmarts('CCC(=O)'), # Keto groups in side chains
        Chem.MolFromSmarts('c1nccc1'), # Uncommon aromatic rings
    ]

    # Check for these unique side chain patterns
    for pattern in unique_modifications:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains unique/modified side chain structures indicating non-proteinogenic nature"

    return False, "Does not have unique/modified side chain features expected in non-proteinogenic amino acids"