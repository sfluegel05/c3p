"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a peptide containing a relatively small number of amino acids (usually <20).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide bonds (-C(=O)NH- patterns)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Count number of amide bonds to estimate chain length
    num_amides = len(amide_matches)
    
    # Each amide bond typically corresponds to one amino acid in the chain (except for the terminal residues)
    num_amino_acids = num_amides + 1  # Assuming a linear chain, last amino acid may not form an additional amide bond
    
    # Threshold for oligopeptides: typically less than 20 amino acids
    if num_amino_acids < 20:
        return True, f"Contains {num_amino_acids} amino acids, which is typical of an oligopeptide"
    else:
        return False, f"Contains {num_amino_acids} amino acids, exceeding typical oligopeptide range"

    return False, "Undefined classification logic"