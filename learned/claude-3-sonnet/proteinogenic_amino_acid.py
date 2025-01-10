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

    # Special case for glycine (achiral)
    glycine_pattern = Chem.MolFromSmarts("[NX3H2][CH2][CX3](=[OX1])[OX2H]")
    if mol.HasSubstructMatch(glycine_pattern):
        # Count heavy atoms (excluding H and isotopes)
        heavy_atoms = set(atom.GetSymbol() for atom in mol.GetAtoms())
        if heavy_atoms.issubset({'C', 'N', 'O'}):
            return True, "Matches glycine structure"

    # Special case for proline
    proline_pattern = Chem.MolFromSmarts("[NX3H0]1[C@H](C(=O)[OH])[CH2][CH2][CH2]1")
    if mol.HasSubstructMatch(proline_pattern):
        return True, "Matches L-proline structure"

    # Special case for pyrrolysine
    pyrrolysine_pattern = Chem.MolFromSmarts("[NX3H2][C@@H](CCCCNC(=O)[C@H]1[C@@H](CC=N1)C)C(=O)[OH]")
    if mol.HasSubstructMatch(pyrrolysine_pattern):
        return True, "Matches L-pyrrolysine structure"

    # Check for basic L-amino acid structure with correct chirality
    l_aa_pattern = Chem.MolFromSmarts("[NX3H2][C@H]([*])[CX3](=[OX1])[OX2H]")
    if not mol.HasSubstructMatch(l_aa_pattern):
        return False, "Does not match L-amino acid structure"

    # Define allowed side chains with strict stereochemistry
    allowed_patterns = [
        # Basic amino acids
        Chem.MolFromSmarts("[NX3H2][C@H](C)[CX3](=[OX1])[OX2H]"),  # Alanine
        Chem.MolFromSmarts("[NX3H2][C@H]([CH](C)C)[CX3](=[OX1])[OX2H]"),  # Valine
        Chem.MolFromSmarts("[NX3H2][C@H](CC(C)C)[CX3](=[OX1])[OX2H]"),  # Leucine
        Chem.MolFromSmarts("[NX3H2][C@H]([C@H](CC)C)[CX3](=[OX1])[OX2H]"),  # Isoleucine
        
        # Sulfur-containing
        Chem.MolFromSmarts("[NX3H2][C@H](CCSC)[CX3](=[OX1])[OX2H]"),  # Methionine
        Chem.MolFromSmarts("[NX3H2][C@H](CS)[CX3](=[OX1])[OX2H]"),  # Cysteine
        
        # Aromatic
        Chem.MolFromSmarts("[NX3H2][C@H](Cc1ccccc1)[CX3](=[OX1])[OX2H]"),  # Phenylalanine
        Chem.MolFromSmarts("[NX3H2][C@H](Cc1ccc(O)cc1)[CX3](=[OX1])[OX2H]"),  # Tyrosine
        Chem.MolFromSmarts("[NX3H2][C@H](Cc1c[nH]cn1)[CX3](=[OX1])[OX2H]"),  # Histidine
        
        # Hydroxy
        Chem.MolFromSmarts("[NX3H2][C@H](CO)[CX3](=[OX1])[OX2H]"),  # Serine
        Chem.MolFromSmarts("[NX3H2][C@H]([C@H](C)O)[CX3](=[OX1])[OX2H]"),  # Threonine
        
        # Acidic and amides
        Chem.MolFromSmarts("[NX3H2][C@H](CC(=O)O)[CX3](=[OX1])[OX2H]"),  # Aspartic acid
        Chem.MolFromSmarts("[NX3H2][C@H](CCC(=O)O)[CX3](=[OX1])[OX2H]"),  # Glutamic acid
        Chem.MolFromSmarts("[NX3H2][C@H](CC(=O)N)[CX3](=[OX1])[OX2H]"),  # Asparagine
        Chem.MolFromSmarts("[NX3H2][C@H](CCC(=O)N)[CX3](=[OX1])[OX2H]"),  # Glutamine
        
        # Basic side chains
        Chem.MolFromSmarts("[NX3H2][C@H](CCCCN)[CX3](=[OX1])[OX2H]"),  # Lysine
        Chem.MolFromSmarts("[NX3H2][C@H](CCCNC(=N)N)[CX3](=[OX1])[OX2H]"),  # Arginine
    ]

    # Check if molecule matches any of the allowed patterns
    matched = False
    for pattern in allowed_patterns:
        if pattern and mol.HasSubstructMatch(pattern):
            matched = True
            break

    if not matched and not (mol.HasSubstructMatch(glycine_pattern) or 
                           mol.HasSubstructMatch(proline_pattern) or 
                           mol.HasSubstructMatch(pyrrolysine_pattern)):
        return False, "Does not match any proteinogenic amino acid pattern"

    # Additional checks for unwanted modifications
    unwanted_patterns = [
        "[CX3](=O)[OX2][CX4]",    # ester
        "[SX4](=O)(=O)",          # sulfone
        "[B,Si,P,I,At]",          # non-standard atoms
        "[NX3]([CX4])([CX4])",    # tertiary amine
    ]
    
    for pattern in unwanted_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Contains non-proteinogenic modifications"

    return True, "Matches proteinogenic amino acid structure"