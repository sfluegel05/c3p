"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determine if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    
    A non-proteinogenic amino acid is defined here as one with both amino and carboxyl groups,
    but with non-standard modifications in side chains or additional groups not seen in
    the 20 standard amino acids.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a non-proteinogenic amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # SMARTS pattern for amino group - broadened to include secondary amines
    amino_group_pattern = Chem.MolFromSmarts('[NX3;H2,H1,H0]')
    # SMARTS pattern for carboxyl group
    carboxyl_group_pattern = Chem.MolFromSmarts('C(=O)[O-,O]')
    
    # Check for amino group and carboxyl group
    has_amino_group = mol.HasSubstructMatch(amino_group_pattern)
    has_carboxyl_group = mol.HasSubstructMatch(carboxyl_group_pattern)
    
    if not (has_amino_group and has_carboxyl_group):
        return False, "Must contain both amino and carboxyl groups"

    # More nuanced detection of unique side chains:
    
    # SMARTS patterns indicative of non-proteinogenic modifications
    modified_side_chain_patterns = [
        # Look for hydroxylated chains, ring alterations, or other unique elements
        Chem.MolFromSmarts('C=C'),    # pi bond in side chain
        Chem.MolFromSmarts('N=C'),    # Amide or imine modifications
        Chem.MolFromSmarts('S'),      # Sulfur in non-cysteine/methionine positions
        Chem.MolFromSmarts('[F,Cl,Br,I]'), # Halogenation as a unique property
        Chem.MolFromSmarts('O=C-N'),  # N-acetyl groups or amidation
        Chem.MolFromSmarts('nc')      # Aromatic nitrogenous bases
    ]

    for pattern in modified_side_chain_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains unique/modified side chain structures indicating non-proteinogenic nature"

    return False, "Does not have unique/modified side chain features expected in non-proteinogenic amino acids"