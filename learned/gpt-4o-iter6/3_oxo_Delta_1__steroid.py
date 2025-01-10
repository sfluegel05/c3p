"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_3_oxo_Delta_1_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    A 3-oxo-Delta(1) steroid is characterized by a steroid structure with
    a double bond between positions 1 and 2 and a ketone group at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # KE each ring can have any stereochemistry, and they need not all be hydrogens).
    # Adjusted SMARTS for steroid backbone with consideration of ring junctions
    steroid_backbone_pattern = "C1C[C@H]2C[C@H](C1)C3C=C(C(C2=O)CC4CCC(=O)C=C3)C4"
    
    # Combined SMARTS pattern for 3-oxo group and Delta 1 double bond within a steroid
    # We need to enforce specific positions for the 3-oxo group and Delta(1) bond
    substructure_pattern = Chem.MolFromSmarts("C1C=C2C3CC[C@@H]4[C@H](C=O)C[C@]4([C@H]3CC=C2C1)")
    
    # Check if the molecule contains the substructure
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(substructure_pattern)):
        return False, "Does not match the 3-oxo-Delta(1) steroid pattern"
    
    # SMARTS pattern Verifying stereochemistry remains flexible for the relevant chiral centers without demanding unnecessary \@/ in SMARTS.
    stereo_pattern = "C1C[C@@H]2C[C@H](C1)C3C=C(/C=O)[C@@]4([C@H]3CCCC4)C2"
    stereo_match = mol.HasSubstructMatch(Chem.MolFromSmarts(stereo_pattern))
    if not stereo_match:
        return False, "Stereochemistry mismatch in the steroid scaffold"

    return True, "Matches the 3-oxo-Delta(1) steroid pattern"