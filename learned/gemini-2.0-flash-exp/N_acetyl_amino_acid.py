"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is an amino acid where the amine group is acetylated
    (CH3-C(=O)-N).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Relaxed SMARTS pattern for N-acetyl amino acid core. The pattern uses a single bond (`~`)
    # to connect to groups to the core structure.
    acetyl_amino_acid_core_pattern = Chem.MolFromSmarts("[CH3][CX3](=[OX1])[NX3]~[CX4]~[CX3](=[OX1])[OX1,OX2-]")

    # Check for the presence of the relaxed core substructure
    if not mol.HasSubstructMatch(acetyl_amino_acid_core_pattern):
        return False, "Molecule does not match the N-acetyl-amino acid core pattern"

    # Check for the number of acetyl groups to be one and only one
    acetyl_pattern = Chem.MolFromSmarts("[CH3][CX3](=[OX1])[NX3]")
    acetyl_matches = mol.GetSubstructMatches(acetyl_pattern)

    if len(acetyl_matches) != 1:
        return False, f"Molecule has {len(acetyl_matches)} acetyl groups, must have only 1"

    # Check for a carboxylic acid group somewhere in the structure.  This check is done separately because
    # the core pattern has already found the required carboxylate, and not necessarily at the end of
    # the chain.
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH,O-]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
       return False, "Molecule does not contain a carboxylic acid group"

    return True, "Molecule matches the N-acetyl-amino acid pattern."