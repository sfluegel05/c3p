"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    Lysophosphatidic acids are characterized by one fatty acid chain, a glycerol backbone, 
    and a terminal phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for glycerol backbone with potential stereochemistry variations
    glycerol_pattern = Chem.MolFromSmarts("OC[C@H](O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for a single ester linkage (only one fatty acid chain)
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected 1 ester linkage, found {len(ester_matches)}"

    # Check for a carbon chain typical of fatty acids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8 or c_count > 24:  # Adjusted range to account for variations
        return False, "Carbon count outside typical range for fatty acids in lysophosphatidic acids"

    # Verify absence of additional glycerol-related esters or chiral centers
    if mol.HasSubstructMatch(Chem.MolFromSmarts("OC[C@H](OC=O)CO")):
        return False, "Extra acyl groups found on glycerol backbone"

    return True, "Structure consistent with lysophosphatidic acid: phosphate group, glycerol backbone, ester linkage for one fatty acid"