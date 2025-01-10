"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    A prostaglandin typically contains a 20-carbon skeleton that includes
    a cyclopentane ring and functional groups such as hydroxyl, ketone, and 
    carboxylic acid attached to the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule resembles a prostaglandin structure, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for prostaglandin core cyclopentane structure with appropriate functionality
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No cyclopentane ring found"

    # Count oxygens and key functional groups - separate patterns might be needed for these
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:  # Prostglandins generally have multiple oxygens
        return False, "Too few oxygen atoms for typical prostaglandin"

    # Check carboxylic acid group presence: -COOH
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing carboxylic acid group"

    # Check for stereochemistry, as prostaglandins are known for this
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 2:
        return False, "Insufficient chirality; prostaglandins typically have chiral centers"

    # Check the chain length to ensure it's around 20 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 22:
        return False, f"Number of carbons is {c_count}, expected around 20 for prostaglandins"

    return True, "Structure contains key features of prostaglandins: cyclopentane ring, carboxyl group, and stereochemistry"