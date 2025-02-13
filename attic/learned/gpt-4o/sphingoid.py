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
    
    # Look for characteristic 2-amino-1,3-diol moiety in sphingoids
    amino_alcohol_pattern = Chem.MolFromSmarts("N[C@H](C)CO")
    if not mol.HasSubstructMatch(amino_alcohol_pattern):
        return False, "No characteristic 2-amino-1,3-diol moiety found"

    # Check for minimal long hydrocarbon chain pattern
    chain_pattern = Chem.MolFromSmarts("CCCCCCCC")  # Minimum length of 8 carbon atoms
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No sufficiently long hydrocarbon chain found"

    # Handle unsaturation as optional with accompanying hydroxyl group requirement
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    hydroxyl_present = any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    if not mol.HasSubstructMatch(double_bond_pattern) and not hydroxyl_present:
        return False, "Lacks either unsaturation or a hydroxyl group, at least one is needed"

    # Verify for presence of stereocenters which are common in sphingoids
    stereo_present = any(atom.HasProp('_CIPCode') for atom in mol.GetAtoms())
    if not stereo_present:
        return False, "No stereochemistry typical of sphingoids found"

    return True, "SMILES corresponds to a recognized sphingoid structure"