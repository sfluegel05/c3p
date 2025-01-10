"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    Attempts to recognize common features in lipopolysaccharides such as lipid A substructures,
    oligosaccharide components, and the potential presence of phosphates.

    Note: Due to structural complexity and diversity, this is a probable, not definitive, classification.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule can be a lipopolysaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Attempt to identify diverse types of fatty acid chains
    # More generic pattern allowing diverse linkages
    lipid_patterns = [
        Chem.MolFromSmarts("C(=O)OC[C;R0][C;R0]"),  # Ester-linked fatty acids
        Chem.MolFromSmarts("C(=O)N[C;R0][C;R0]")   # Amide-linked fatty acids
        # Further patterns could be added here
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in lipid_patterns):
        return (False, "No lipid chain motifs detected")

    # Enhanced recognition of saccharide units, include more complex variations
    saccharide_pattern = Chem.MolFromSmarts("[C@@H]1(O)[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O1)CO")  # Mono-/oligosaccharide pattern
    if not mol.HasSubstructMatch(saccharide_pattern):
        return (False, "No polysaccharide-like structures detected")

    # Detect phosphate motifs but less strictly, allowing for various representations
    phosphate_pattern = Chem.MolFromSmarts("[O-]P(=O)(O)O")  # Common phosphate representation
    if not mol.HasSubstructMatch(phosphate_pattern):
        return (False, "No phosphate groups identified; common but not universal in LPS")

    mol_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 600:
        return (False, "Molecular weight too low to likely be an LPS")

    return (True, "Contains structures consistent with lipopolysaccharides: lipids, polysaccharides, and potential phosphates")