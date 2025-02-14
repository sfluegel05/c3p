"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdChemReactions

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for L-alpha-amino acid pattern
    # Carbon with a chiral center having attached NH2 and COOH, and verify L-configuration
    # Also require this not to be part of a peptide chain
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[C@@H](N)[C;!$([N]C=O)](=O)O")
    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "Does not match the L-alpha-amino acid structure"

    # Verify stereochemistry for L-configuration at correct center
    # The COO and NH2 attached chiral center must be (S) stereochemistry in most standards
    chiral_matches = mol.GetSubstructMatches(alpha_amino_acid_pattern)
    for match in chiral_matches:
        atom = mol.GetAtomWithIdx(match[0])
        if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
            continue
        
        # Verify if the pattern occurs isolated, not forming peptide bonds
        neighbors = atom.GetNeighbors()
        num_carboxylic_acid = sum(1 for n in neighbors if n.GetSymbol() == 'O' and n.GetDoubleProp("_CIPRank") > atom.GetDoubleProp("_CIPRank"))
        num_amino = sum(1 for n in neighbors if n.GetSymbol() == 'N')
        
        if num_carboxylic_acid == 2 and num_amino == 1:
            return True, "Contains L-alpha-amino acid structure"

    return False, "Does not match the isolated L-alpha-amino acid structure without peptide linkage"