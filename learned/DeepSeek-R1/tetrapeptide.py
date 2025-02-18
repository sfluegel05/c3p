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
    
    # Improved patterns accounting for linear (3 bonds) and cyclic (4 bonds) structures
    # Linear pattern allows modified termini (e.g. acetylated N-term or amidated C-term)
    linear_pattern = Chem.MolFromSmarts(
        "[$([NH2]-[CH]),$(C(=O)-[NH]-[CH])]"  # N-term
        "[CH](-[!#1])-C(=O)-[NH]-[CH](-[!#1])-C(=O)-[NH]-[CH](-[!#1])-C(=O)-[NH]-[CH](-[!#1])"  # Core 3 bonds
        "-C(=O)-[$([OH]),$([NH2]),$([O-]),$([NH]-[!#1])]"  # C-term
    )
    
    # Cyclic pattern requires 4 peptide bonds in a ring
    cyclic_pattern = Chem.MolFromSmarts(
        "C(=O)-[NH]-[CH](-[!#1])-C(=O)-[NH]-[CH](-[!#1])-C(=O)-[NH]-[CH](-[!#1])-C(=O)-[NH]-[CH](-[!#1])-C(=O)"
    )
    
    # Check for matches
    has_linear = mol.HasSubstructMatch(linear_pattern)
    has_cyclic = mol.HasSubstructMatch(cyclic_pattern)
    
    if has_linear or has_cyclic:
        return True, "Contains four amino acid residues connected by peptide bonds"
    
    # Fallback: Verify peptide bond count and alpha carbons
    peptide_bonds = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=[OX1])-[NX3]")))
    alpha_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CH](-[!#1])(-C(=O)-[NX3])")))
    
    # Linear requires 3 bonds, cyclic 4. Alpha carbons should be >=4 (modified residues may have more)
    if (peptide_bonds == 3 and alpha_carbons >=4) or (peptide_bonds ==4 and alpha_carbons ==4):
        return True, f"Found {peptide_bonds} peptide bonds and {alpha_carbons} alpha carbons"
    
    return False, "Does not match tetrapeptide structural requirements"