"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: CHEBI:27480 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol has the general formula H-[CH2C(Me)=CHCH2]nOH, where the carbon skeleton
    is composed of more than one isoprene unit.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for isoprene unit pattern (C=C-C-C=C)
    isoprene_pattern = Chem.MolFromSmarts("[CH2]=[C]([CH3])[CH2][CH2][CH2]=[C]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    # Check for at least two isoprene units
    if len(isoprene_matches) < 2:
        return False, "Less than two isoprene units found"
    
    # Check for terminal hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_match = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_match:
        return False, "No terminal hydroxyl group found"
    
    # Check for a hydrocarbon skeleton with possible additional functional groups
    heteroatom_pattern = Chem.MolFromSmarts("[!#6;!#1;!#8]")  # Non-carbon, non-hydrogen, non-oxygen atoms
    heteroatom_matches = mol.GetSubstructMatches(heteroatom_pattern)
    if len(heteroatom_matches) > 0:
        return False, "Found heteroatoms other than oxygen"
    
    # Check for a minimum carbon chain length (e.g., 10 carbon atoms)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 10:
        return False, "Carbon chain too short for a polyprenol"
    
    return True, "Molecule has the general formula H-[CH2C(Me)=CHCH2]nOH with a hydrocarbon skeleton"

# Example usage
smiles = "OC(CCC=C(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(=CCCC(=CCCC(=CCCC(=CCO)C)C)C)C)C)C)C)C"
print(is_polyprenol(smiles))