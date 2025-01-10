"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is characterized as an aliphatic monocarboxylic acid with a chain of 4 to 28 carbons, 
    which may be saturated or unsaturated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Get the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Fatty acids typically have between 4 and 28 carbon atoms
    if c_count < 4 or c_count > 28:
        return False, f"Carbon chain length {c_count} not in [4, 28]"

    # Ensure the molecule is essentially aliphatic
    has_ring = mol.GetRingInfo().NumRings() > 0
    # Fatty acids should not have rings by definition
    if has_ring:
        return False, "Contains ring structure(s)"
    
    # Check the acyclic nature of carbon chain
    main_chain = max(rdMolDescriptors.CalcMolecularGraphDescriptors(mol, "NumAtoms") or [0], key=len)
    if len(main_chain) < c_count - 1:
        return False, "Main chain is too short, likely the molecule is branched"

    return True, "Contains carboxylic acid group with an aliphatic chain between 4 and 28 carbons"