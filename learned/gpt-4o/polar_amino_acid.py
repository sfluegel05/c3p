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

    # Check for alpha-amino acid structure
    amino_acid_pattern = Chem.MolFromSmarts("[N][C@@H,C@H](C(O)=O)")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Not an alpha-amino acid"

    # SMARTS patterns for hydrogen-bonding functional groups in side chains
    polar_side_chain_patterns = [
        Chem.MolFromSmarts("[CX3](=O)[NX3]"),  # Amide group - Asparagine, Glutamine
        Chem.MolFromSmarts("[OX2H]"),          # Hydroxyl group - Serine, Threonine, Tyrosine 
        Chem.MolFromSmarts("[NX3][CX3](=[OX1])"),  # Guanidinium group - Arginine
        Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]"),# Carboxyl group - Aspartic acid, Glutamic acid
        Chem.MolFromSmarts("[#16]"),           # Thiol group - Cysteine
    ]

    # Check for polar side chain patterns
    for pattern in polar_side_chain_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Polar side chain functional group found"
    
    return False, "No polar side chain functional groups found"

# Test the function with example polar amino acids
example_smiles = [
    "N[C@H](CS)C(O)=O",  # D-cysteine
    "NC(Cc1c[nH]c2ccccc12)C(O)=O",  # tryptophan
    "NC(Cc1ccc(O)cc1)C(O)=O",  # tyrosine
]

for smi in example_smiles:
    result, reason = is_polar_amino_acid(smi)
    print(f"SMILES: {smi}, Polar: {result}, Reason: {reason}")