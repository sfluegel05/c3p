"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids include sphinganine, its homologs and stereoisomers, and hydroxyl or unsaturated derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the characteristic amino alcohol moiety in sphingoids
    amino_alcohol_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#6](N)-[#6](O)-[#6]")  # Simple pattern for amino alcohol
    if not mol.HasSubstructMatch(amino_alcohol_pattern):
        return False, "No characteristic amino alcohol moiety found"
    
    # Check for long hydrocarbon chain, possibly branched
    chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCC")  # Long hydrocarbon chain as an example
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long hydrocarbon chains found"

    # Optional check for backbone stereochemistry typically seen in sphingoids
    stereo_centers = [atom.GetChiralTag() for atom in mol.GetAtoms()]
    has_stereo = Chem.CHI_TetrahedralIsomer.Any in stereo_centers
    if not has_stereo:
        return False, "No stereochemistry typical of sphingoids found"
    
    return True, "SMILES corresponds to a recognized sphingoid structure"