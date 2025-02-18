"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
"""
Classifies: CHEBI_83820 non-proteinogenic amino acid
"""
from rdkit import Chem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    Non-proteinogenic amino acids are those not naturally encoded in the genetic code.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check alpha-amino acid structure: amine adjacent to carboxyl group
    alpha_amino_pattern = Chem.MolFromSmarts('[#7]C-C(=O)O')
    if not mol.HasSubstructMatch(alpha_amino_pattern):
        return False, "Not an alpha-amino acid"
    
    # Generate canonical SMILES with stereochemistry
    input_canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    
    # Predefined standard L-amino acids (canonical SMILES with stereochemistry)
    standard_amino_smiles = [
        'C(C(=O)O)N',  # Glycine
        'C[C@H](N)C(=O)O',  # L-Alanine
        'CC(C)[C@H](N)C(=O)O',  # L-Valine
        'CC(C)CC[C@H](N)C(=O)O',  # L-Leucine
        'CCC[C@H](C)[C@H](N)C(=O)O',  # L-Isoleucine
        'N1[C@H](CCC1)C(=O)O',  # L-Proline
        'C([C@H](N)C(=O)O)O',  # L-Serine
        'C[C@H](O)[C@H](N)C(=O)O',  # L-Threonine
        'C([C@H](N)C(=O)O)S',  # L-Cysteine
        'CSCC[C@H](N)C(=O)O',  # L-Methionine
        'C([C@H](N)C(=O)O)C(=O)N',  # L-Asparagine
        'C(CC(=O)N)[C@H](N)C(=O)O',  # L-Glutamine
        'C(CCCCN)[C@H](N)C(=O)O',  # L-Lysine
        'C(CCCNC(=N)N)[C@H](N)C(=O)O',  # L-Arginine
        'C1=C(NC=N1)C[C@H](N)C(=O)O',  # L-Histidine
        'C([C@H](N)C(=O)O)C(=O)O',  # L-Aspartic acid
        'C(CC(=O)O)[C@H](N)C(=O)O',  # L-Glutamic acid
        'C1=CC=CC=C1C[C@H](N)C(=O)O',  # L-Phenylalanine
        'C1=CC(=CC=C1C[C@H](N)C(=O)O)O',  # L-Tyrosine
        'C1=CC=C2C(=C1)C(=CN2)C[C@H](N)C(=O)O'  # L-Tryptophan
    ]
    
    # Generate canonical SMILES for each standard amino acid
    standard_canonical = set()
    for s in standard_amino_smiles:
        m = Chem.MolFromSmiles(s)
        if m is not None:
            canonical = Chem.MolToSmiles(m, isomericSmiles=True)
            standard_canonical.add(canonical)
    
    if input_canonical in standard_canonical:
        return False, "Standard proteinogenic amino acid"
    else:
        return True, "Non-standard alpha-amino acid structure"