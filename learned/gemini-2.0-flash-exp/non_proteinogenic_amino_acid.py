"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    A non-proteinogenic amino acid is an amino acid not naturally encoded in the genetic code.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Basic amino acid structure check: N-C-C(=O)-O and a sidechain
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3][CX4]([H])[CX3](=[OX1])[OX2]")
    if not mol.HasSubstructMatch(alpha_carbon_pattern):
        return False, "Missing basic amino acid structure"

    # 2. List of SMILES for the 20 standard proteinogenic amino acids (L-forms, and also D for Asp and Glu)
    canonical_amino_acids_smiles = [
        "N[C@@H](C)C(O)=O",  # Alanine
        "N[C@@H](CC(C)=O)C(O)=O", # Asparagine
        "N[C@@H](CC(=O)O)C(O)=O", # Aspartic acid
        "N[C@@H](Cc1c[nH]cn1)C(O)=O",  # Histidine
        "N[C@@H](Cc1ccccc1)C(O)=O",  # Phenylalanine
        "N[C@@H](CS)C(O)=O", # Cysteine
        "N[C@@H](CC(=O)N)C(O)=O", # Glutamine
        "N[C@@H](CCC(=O)O)C(O)=O", # Glutamic acid
        "N[C@@H](C(C)C)C(O)=O", # Valine
        "N[C@@H](CCCNC(=N)N)C(O)=O", # Arginine
        "N[C@@H](CO)C(O)=O", # Serine
        "N[C@@H]([C@@H](O)C)C(O)=O",  # Threonine
        "N[C@@H](CCSC)C(O)=O",  # Methionine
        "N[C@@H](C(C)CC)C(O)=O", # Leucine
        "N[C@@H](C(C)(C)C)C(O)=O",  # Isoleucine
        "N[C@@H](Cc1c[nH]c2ccccc12)C(O)=O", # Tryptophan
        "N[C@@H](Cc1ccc(O)cc1)C(O)=O",  # Tyrosine
        "N1C[C@@H]2[C@H](C[C@@H]1C(O)=O)NC2",  # Proline
        "N[C@H](CC(O)=O)C(O)=O", # D-aspartic acid
        "N[C@H](CCC(O)=O)C(O)=O" # D-glutamic acid
    ]
    
    # Canonical check. Convert to canonical SMILES for better comparison
    canonical_smiles = [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in canonical_amino_acids_smiles]
    input_smiles = Chem.MolToSmiles(mol)

    if input_smiles in canonical_smiles:
        return False, "It is a standard proteinogenic amino acid."
    
    return True, "It is a non-proteinogenic amino acid."