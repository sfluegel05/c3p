"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
"""
Classifies: CHEBI_83820 non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import Mol, MolFromSmiles, MolFromSmarts
from rdkit.Chem.rdmolops import GetFormalCharge

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid.
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for alpha-amino acid structure: NH2-CH(R)-COOH
    alpha_amino_pattern = MolFromSmarts('[NH2,NX3H,NX4H2]-[CX4H]([!O])-[CX3](=O)[OH]')
    if not mol.HasSubstructMatch(alpha_amino_pattern):
        return False, "Not an alpha-amino acid"
    
    # Check for peptide bonds (amide linkages between amino acids)
    peptide_pattern = MolFromSmarts('[CX3](=O)-[NX3H]-[CX4]')
    if mol.HasSubstructMatch(peptide_pattern):
        return False, "Contains peptide bond"
    
    # Standard amino acid core structures (exact matches)
    standard_amino_smarts = [
        # Glycine
        '[NH2]C([H])(C(=O)O)',
        # L-Alanine
        'C[C@H]([NH2])C(=O)O',
        # L-Valine
        'CC(C)[C@H]([NH2])C(=O)O',
        # L-Leucine
        'CC(C)CC[C@H]([NH2])C(=O)O',
        # L-Isoleucine
        'CCC[C@H](C)[C@H]([NH2])C(=O)O',
        # L-Proline (cyclic)
        'N1C([CH2][CH2]C1)C(=O)O',
        # L-Serine
        'C([C@H]([NH2])C(=O)O)O',
        # L-Threonine
        'C[C@H](O)[C@H]([NH2])C(=O)O',
        # L-Cysteine
        'C([C@H]([NH2])C(=O)O)S',
        # L-Methionine
        'CSCC[C@H]([NH2])C(=O)O',
        # L-Asparagine
        'C([C@H]([NH2])C(=O)O)C(=O)N',
        # L-Glutamine
        'C(CC(=O)N)[C@H]([NH2])C(=O)O',
        # L-Lysine
        'C(CCCCN)[C@H]([NH2])C(=O)O',
        # L-Arginine
        'C(CCCNC(=N)N)[C@H]([NH2])C(=O)O',
        # L-Histidine
        'C1=C(NC=N1)C[C@H]([NH2])C(=O)O',
        # L-Aspartic acid
        'C([C@H]([NH2])C(=O)O)C(=O)O',
        # L-Glutamic acid
        'C(CC(=O)O)[C@H]([NH2])C(=O)O',
        # L-Phenylalanine
        'C1=CC=CC=C1C[C@H]([NH2])C(=O)O',
        # L-Tyrosine
        'C1=CC(=CC=C1C[C@H]([NH2])C(=O)O)O',
        # L-Tryptophan
        'C1=CC=C2C(=C1)C(=CN2)C[C@H]([NH2])C(=O)O'
    ]
    
    # Check if the entire molecule matches any standard amino acid
    for smarts in standard_amino_smarts:
        pattern = MolFromSmarts(smarts)
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return False, "Matches standard amino acid"
    
    # Check for zwitterionic forms of standard amino acids
    zwitterion_pattern = MolFromSmarts('[NH3+]C([H])(C(=O)[O-])')
    if mol.HasSubstructMatch(zwitterion_pattern):
        return False, "Standard amino acid zwitterion"
    
    # If passed all checks, likely non-proteinogenic
    return True, "Non-standard amino acid structure"