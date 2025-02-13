"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride contains a glycerol backbone with one fatty acid chain
    attached via an ester bond at the 1-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern for glycerol backbone: OH groups at 1, 2, 3 positions
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Pattern for 1-position ester linkage: ester linkage with carbon chain
    ester_pattern = Chem.MolFromSmarts("C(=O)OC(CO)CO")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found at expected position"
    
    # Check for only one ester linkage
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)O"))
    if len(ester_matches) != 1:
        return False, f"Expected 1 ester linkage, found {len(ester_matches)}"
    
    # Consideration of chirality, technically 1-monoglycerides can be racemic
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 1:
        return False, f"Expected at least one chiral center, found {len(chiral_centers)}"
    
    return True, "Contains glycerol backbone with one fatty acid chain attached via ester bond at position 1"