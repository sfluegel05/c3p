"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid has a side chain capable of forming one or more hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the basic alpha-amino acid backbone
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[N;!R][C@H,C@@H](C(O)=O)")
    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "Not an alpha-amino acid"
    
    # SMARTS patterns for hydrogen-bonding functional groups in side chains
    polar_side_chain_patterns = [
        Chem.MolFromSmarts("NC(=O)[#6]"),      # Amide side chain - Asparagine or Glutamine
        Chem.MolFromSmarts("[OH][CH2,C]"),     # Alcohol side chain - Serine or Threonine
        Chem.MolFromSmarts("[nH]"),            # Indole or imidazole group - Tryptophan or Histidine
        Chem.MolFromSmarts("[NX3H2]"),         # Ammonium primary amine - Lysine
        Chem.MolFromSmarts("NC(N)=N"),         # Guanidinium - Arginine
        Chem.MolFromSmarts("[SH]"),            # Thiol group - Cysteine
        Chem.MolFromSmarts("ccc(O)"),          # Phenol group - Tyrosine
        Chem.MolFromSmarts("C(=O)[OX1H]")      # Carboxyl group - Glutamic acid, Aspartic acid
    ]

    # Check for polar side chain patterns
    for pattern in polar_side_chain_patterns:
        if pattern and mol.HasSubstructMatch(pattern):  # Ensure the pattern is valid
            return True, "Polar side chain functional group found"
    
    return False, "No polar side chain functional groups found"

# Sample Test
example_smiles = [
    "N[C@H](CS)C(O)=O",  # D-cysteine
    "NC(Cc1c[nH]c2ccccc12)C(O)=O",  # tryptophan
    "NC(Cc1ccc(O)cc1)C(O)=O",  # tyrosine
]

for smi in example_smiles:
    result, reason = is_polar_amino_acid(smi)
    print(f"SMILES: {smi}, Polar: {result}, Reason: {reason}")