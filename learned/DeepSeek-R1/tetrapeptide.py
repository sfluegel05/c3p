"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide contains four amino acid residues connected by peptide bonds in a linear or cyclic chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern to match peptide backbone with four residues
    # Accounts for linear (3 bonds) and cyclic (4 bonds) structures
    # The pattern matches four amino groups connected via alpha carbons and carbonyls
    tetrapeptide_pattern = Chem.MolFromSmarts(
        "[NH]-[CH](-[!#1])(-[!#1])-C(=O)-[NH]-[CH](-[!#1])(-[!#1])-C(=O)-[NH]-[CH](-[!#1])(-[!#1])-C(=O)-[NH]-[CH](-[!#1])(-[!#1])-C(=O)"
    )
    cyclic_pattern = Chem.MolFromSmarts(
        "C(=O)-[NH]-[CH](-[!#1])(-[!#1])-C(=O)-[NH]-[CH](-[!#1])(-[!#1])-C(=O)-[NH]-[CH](-[!#1])(-[!#1])-C(=O)-[NH]-[CH](-[!#1])(-[!#1])-C(=O)"
    )
    
    # Check for linear tetrapeptide (3 peptide bonds) or cyclic (4 bonds)
    has_linear = mol.HasSubstructMatch(tetrapeptide_pattern)
    has_cyclic = mol.HasSubstructMatch(cyclic_pattern)
    
    if has_linear or has_cyclic:
        return True, "Contains four amino acid residues connected by peptide bonds"
    
    # Additional check for modified termini (e.g., acetylated N-terminus or amidated C-terminus)
    modified_linear_pattern = Chem.MolFromSmarts(
        "[C;H0](=O)-[NH]-[CH](-[!#1])(-[!#1])-C(=O)-[NH]-[CH](-[!#1])(-[!#1])-C(=O)-[NH]-[CH](-[!#1])(-[!#1])-C(=O)-[NH]-[CH](-[!#1])(-[!#1])-C(=O)-[OH0]"
    )
    if mol.HasSubstructMatch(modified_linear_pattern):
        return True, "Modified tetrapeptide backbone detected"
    
    # Count peptide bonds as fallback (may catch some edge cases)
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])-[NX3]")
    num_peptide_bonds = len(mol.GetSubstructMatches(peptide_bond_pattern))
    
    if num_peptide_bonds >= 3 and num_peptide_bonds <=4:
        # Verify there are exactly four amino acid residues
        amino_groups = rdMolDescriptors.CalcNumAmideBonds(mol)
        if amino_groups >=4:
            return True, f"Found {num_peptide_bonds} peptide bonds and four amino groups"
    
    return False, "Does not match tetrapeptide structural requirements"