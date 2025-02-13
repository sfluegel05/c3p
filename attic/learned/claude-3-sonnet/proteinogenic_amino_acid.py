"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: proteinogenic amino acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for peptide bonds
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4][NX3][CX3](=[OX1])")
    if mol.HasSubstructMatch(peptide_pattern):
        return False, "Contains peptide bonds"

    # Check for N-modifications (except for proline)
    n_modification_pattern = Chem.MolFromSmarts("[CX4][NX3]([!H])[CX3](=O)[OX2H]")
    if mol.HasSubstructMatch(n_modification_pattern):
        return False, "Contains N-modifications"

    # Special case for glycine (achiral)
    glycine_pattern = Chem.MolFromSmarts("[NX3H2][CH2][CX3](=[OX1])[OX2H]")
    if mol.HasSubstructMatch(glycine_pattern):
        # Check if it's a simple glycine (allowing isotopes)
        atom_symbols = set(a.GetSymbol() for a in mol.GetAtoms())
        if atom_symbols.issubset({'C', 'H', 'N', 'O'}):
            return True, "Matches glycine structure"

    # Check for basic L-amino acid structure with correct chirality
    l_aa_pattern = Chem.MolFromSmarts("[NX3H2][C@H][CX3](=[OX1])[OX2H]")
    l_pro_pattern = Chem.MolFromSmarts("[NX3H0]1[C@H](C(=O)[OH])[CH2][CH2][CH2]1")
    
    is_l_aa = mol.HasSubstructMatch(l_aa_pattern)
    is_proline = mol.HasSubstructMatch(l_pro_pattern)
    
    if not (is_l_aa or is_proline):
        return False, "Does not match L-amino acid structure"

    # Define allowed side chains with stereochemistry where relevant
    side_chains = {
        'alanine': '[C@H](C)(N)C(=O)O',
        'valine': '[C@H]([CH](C)C)(N)C(=O)O',
        'leucine': '[C@H](CC(C)C)(N)C(=O)O',
        'isoleucine': '[C@H]([C@H](CC)C)(N)C(=O)O',
        'proline': 'C1CC[C@@H](N1)C(=O)O',
        'methionine': '[C@H](CCSC)(N)C(=O)O',
        'phenylalanine': '[C@H](Cc1ccccc1)(N)C(=O)O',
        'serine': '[C@H](CO)(N)C(=O)O',
        'threonine': '[C@H]([C@H](C)O)(N)C(=O)O',
        'cysteine': '[C@H](CS)(N)C(=O)O',
        'tyrosine': '[C@H](Cc1ccc(O)cc1)(N)C(=O)O',
        'asparagine': '[C@H](CC(=O)N)(N)C(=O)O',
        'glutamine': '[C@H](CCC(=O)N)(N)C(=O)O',
        'aspartic acid': '[C@H](CC(=O)O)(N)C(=O)O',
        'glutamic acid': '[C@H](CCC(=O)O)(N)C(=O)O',
        'lysine': '[C@H](CCCCN)(N)C(=O)O',
        'arginine': '[C@H](CCCNC(=N)N)(N)C(=O)O',
        'histidine': '[C@H](Cc1c[nH]cn1)(N)C(=O)O',
        'selenocysteine': '[C@H](C[SeH])(N)C(=O)O',
        'pyrrolysine': '[C@H](CCCCNC(=O)[C@H]1CC=NC1C)(N)C(=O)O'
    }

    # Check if molecule matches any proteinogenic amino acid core structure
    matched = False
    for aa, pattern in side_chains.items():
        # Convert pattern to mol with explicit H to match isotope-labeled compounds
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol and mol.HasSubstructMatch(pattern_mol):
            matched = True
            break

    if not matched and not mol.HasSubstructMatch(glycine_pattern):
        return False, "Side chain doesn't match any proteinogenic amino acid"

    # Additional check for unwanted groups
    unwanted_groups = [
        "[NX3]([CX4])[CX3](=O)",  # peptide bond
        "[CX3](=O)[OX2][CX4]",    # ester
        "[CX4][SX4](=O)(=O)",     # sulfone
        "[B,Si,P,I,At]"           # non-standard atoms
    ]
    
    for pattern in unwanted_groups:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Contains non-proteinogenic modifications"

    return True, "Matches proteinogenic amino acid structure"