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

    # Pattern for isoprene unit (repeated in polyprenols): [C=C(C)CCC]
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)CCC")
    isoprene_matches = len(mol.GetSubstructMatches(isoprene_pattern))
    if isoprene_matches < 3:
        # At least 3 isoprene units expected for polyprenyl chains
        return False, f"Found {isoprene_matches} isoprene units, need at least 3 for a polyprenol chain"

    # Pattern for phosphate group: [OP(=O)(O)O]
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 1:
        return False, "No phosphate group found"

    # Check for connection between polyprenol chain and phosphate group
    if not mol.HasSubstructMatch(Chem.MolFromSmarts("CCC[OH]P(=O)(O)O")):
        # This is a strong pattern check combining a polyprenol terminal with phosphate linkage
        return False, "Polyprenol chain not appropriately connected to phosphate group"

    return True, "Molecule is a polyprenol phosphate"

# Example usage:
# smiles = "CC(C)=CCC\C(C)=C\CC\C(C)=C\CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/COP(O)(=O)OP(O)(O)=O"
# print(is_polyprenol_phosphate(smiles))