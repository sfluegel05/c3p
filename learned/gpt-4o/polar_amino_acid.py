"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Check for basic alpha-amino acid backbone
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[N;!R][C@H,C@@H](C(O)=O)")
    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "Not an alpha-amino acid"
    
    # SMARTS patterns for hydrogen-bonding functional groups in side chains
    # Eliminate patterns that are generally found on the alpha backbone
    polar_side_chain_patterns = [
        Chem.MolFromSmarts("[$([NX3][CX3](=O)[#6])!H0]"),  # Amide side chain - Asparagine, Glutamine
        Chem.MolFromSmarts("[$([OH][^C](C)[O])!H0]"),       # Alcohol side group - Serine, Threonine, Tyrosine
        Chem.MolFromSmarts("[$([NX3H2][CX3]=N[CX3]=N)!H0]"),# Guanidinium - Arginine
        Chem.MolFromSmarts("[$([CX3](=O)[OX2H1]!C){}]"),    # Unusual carboxylate - ignore if in backbone
        Chem.MolFromSmarts("[$([S][#6])!H0]"),              # Thiol group - Cysteine
    ]

    non_backbone_oxygen = Chem.MolFromSmarts("[CX3](=O)[OX2H1]!([$([C][OH]=O),$([C]CC[N]!N)])") # Exclude backbone

    # Check for polar side chain patterns specifically
    for pattern in polar_side_chain_patterns:
        if mol.HasSubstructMatch(pattern):
            if mol.HasSubstructMatch(non_backbone_oxygen):
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