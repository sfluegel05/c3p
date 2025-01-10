"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is a glycerol moiety with a single acyl group attached at
    the primary position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1-monoglyceride, False otherwise
        str: Explanation of the result
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone (we recognize primary alcohols and 1 secondary alcohol)
    glycerol_backbone = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_backbone):
        return False, "No glycerol backbone found"

    # Elicit an expanded check for ester linkage that may allow for common stereochemical variations
    ester_primary_position = Chem.MolFromSmarts("C(O)C(=O)[*]")
    ester_matches = mol.GetSubstructMatches(ester_primary_position)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1 at the primary position"

    # Chiral check - relaxed to count stereocenters but not demand specific stereo SMILES unless specifically sn-stereo is indicated 
    chiral_center_check = any(atom.GetChiralTag() for atom in mol.GetAtoms())
    if not chiral_center_check:
        return True, "Contains glycerol backbone with one acyl linkage but no stereochemistry is detected uniquely (assuming achiral or racemic mix)"

    # Validated based on stereo if required by the context (optional, as some could be achiral)
    sn_glycerol = Chem.MolFromSmarts("[C@H](O)C(O)[CH2]")
    if any(mol.HasSubstructMatch(sn_glycerol)):
        return True, "Contains sn-glycerol backbone with a single acyl group esterified at a primary position"
    
    return False, "Failed to precisely conform to a recognized 1-monoglyceride structure"