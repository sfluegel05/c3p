"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
# Import required libraries from RDKit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    Monoterpenoid indole alkaloids are complex structures derived from tryptophan and specialized 
    diisoprenoid precursors, with characteristic indole and expanded heterocyclic rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded indole-like core structure patterns including fusion possibilities
    indole_patterns = [
        Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3ccccc23"),  # Standard indole
        Chem.MolFromSmarts("c1cc2[nH]c3ccccc3c2c1"),  # Fused indole
        Chem.MolFromSmarts("c1cc2cnc(c2cc1)N"),  # Indolinamine
        Chem.MolFromSmarts("c1cc2[nH]cnc2c1"),  # Pyrroloindoles
        Chem.MolFromSmarts("c2cc1[nH]c3ccccc3nc1cc2")  # Indicative of complex indole systems
    ]
    if not any(mol.HasSubstructMatch(pat) for pat in indole_patterns):
        return False, "No recognizable indole-like structure"

    # Count nitrogen atoms - ensure sufficient number for alkaloids
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:  # Reflecting the potential complexity and additional nitrogen heterocycles
        return False, "Insufficient nitrogen atoms for typical alkaloid complexity"

    # Total ring count and aromatic complexity
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_rings < 3 or num_aromatic_rings < 2:  # Expecting more complex and aromatic systems
        return False, f"Not enough total or aromatic rings: {num_rings} total, {num_aromatic_rings} aromatic"

    # Check for monoterpenoid-like feature
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)C")
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "Missing diisoprenoid units typical of monoterpenoids"

    # Check for other functional group patterns specific to monoterpenoid indole alkaloids
    # (methoxy enhancements have been implied but can involve other groups—in-depth specificity could be added here)

    return True, "Recognizable complex indole structure with monoterpenoid features and sufficient nitrogen"