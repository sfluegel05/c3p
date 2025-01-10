"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a given SMILES string corresponds to a polyprenol phosphate.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for isoprene units, generalized pattern to capture stereochemistry and variations
    # Generically detect: [C=C(C)CC]
    isoprene_pattern_general = Chem.MolFromSmarts("C=C(C)CC")
    isoprene_matches = len(mol.GetSubstructMatches(isoprene_pattern_general))
    if isoprene_matches < 3:
        return False, f"Found {isoprene_matches} isoprene units, need at least 3 for a polyprenol chain"
    
    # Identify phosphate group, adjusted pattern: [OP(=O)(O)]
    phosphate_pattern_general = Chem.MolFromSmarts("OP(=O)(O)")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern_general))
    if phosphate_matches < 1:
        return False, "No phosphate group found"

    # Check for connection between polyprenol chain and phosphate group (terminal C to P linkage)
    # Relax strict terminal checking, check extension from one end involves O-P linkage
    connection_pattern = Chem.MolFromSmarts("C-O-P(=O)(O)C") 
    if not mol.HasSubstructMatch(connection_pattern):
        return False, "Polyprenol chain not appropriately connected to phosphate group"

    return True, "Molecule is a polyprenol phosphate"

# Example usage:
# smiles = "CC(C)=CCC\C(C)=C\CC\C(C)=C\CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/COP(O)(=O)OP(O)(O)=O"
# print(is_polyprenol_phosphate(smiles))