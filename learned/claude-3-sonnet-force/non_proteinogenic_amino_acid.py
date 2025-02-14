"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
"""
Classifies: CHEBI:33673 non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

PROTEINOGENIC_SMILES = ['NC(C(=O)O)C(C)C',  # Alanine
                         'NC(C(=O)O)CC(C)C',  # Valine
                         'NC(C(=O)O)CC(=O)O',  # Aspartic acid
                         'NC(C(=O)O)C(C)(C)C',  # Leucine
                         'NC(C(=O)O)CC1=CNC=N1',  # Histidine
                         'NC(C(=O)O)CS',  # Cysteine
                         'NC(C(=O)O)C(O)C(=O)O',  # Glutamic acid
                         'NC(C(=O)O)C1=CN=CN1',  # Tryptophan
                         'NC(C(=O)O)CO',  # Serine
                         'NC(C(=O)O)CC1=CC=C(O)C=C1',  # Tyrosine
                         'NC(C(=O)O)C(N)=N',  # Arginine
                         'NC(C(=O)O)CC1=CNc2ccccc12',  # Tryptophan
                         'NC(C(=O)O)CCSC(N)=N',  # Methionine
                         'NC(C(=O)O)CC(C)C(=O)N',  # Proline
                         'NC(C(=O)O)CC(=O)NC(N)=N',  # Citrulline
                         'NC(C(=O)O)CC1=CNC2=C1C=CC=C2',  # Histidine
                         'NC(C(=O)O)CC(N)=O',  # Asparagine
                         'NC(C(=O)O)C(C)(C)O',  # Threonine
                         'NC(C)C(=O)O',  # Glycine
                         'NC(C(=O)O)CC1=CC=CC=C1',  # Phenylalanine
                         'NC(C(=O)O)CCCNC(N)=N',  # Lysine
                         'NC(C(=O)O)CCC(C)C',  # Isoleucine
                         'NC(C(=O)NCC(N)=O)C(C)C',  # Glutamine
                         ]


def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    A non-proteinogenic amino acid is any amino acid that is not naturally encoded in the genetic code.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carboxyl and amino groups
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    amino_pattern = Chem.MolFromSmarts("N[C;H3]")
    n_carboxyl = len(mol.GetSubstructMatches(carboxyl_pattern))
    n_amino = len(mol.GetSubstructMatches(amino_pattern))

    # Amino acid must have one carboxyl and one amino group
    if n_carboxyl != 1 or n_amino != 1:
        return False, "Does not contain one carboxyl and one amino group"

    # Check if it matches a proteinogenic amino acid
    is_proteinogenic = False
    for smiles in PROTEINOGENIC_SMILES:
        proteinogenic_mol = Chem.MolFromSmiles(smiles)
        if mol.HasSubstructMatch(proteinogenic_mol):
            is_proteinogenic = True
            break

    if is_proteinogenic:
        return False, "Matches a proteinogenic amino acid"
    else:
        return True, "Contains carboxyl and amino groups but does not match any proteinogenic amino acids"