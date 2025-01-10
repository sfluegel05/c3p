"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
from rdkit import Chem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a heterocyclic fatty acid containing an epoxide ring as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for epoxide ring pattern, which is a three-membered ring with an oxygen
    epoxide_pattern = Chem.MolFromSmarts("[C]-[O]-[C]")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide ring found"

    # Check for a long carbon chain characteristic of fatty acids
    # This includes the constraint for length and branched structure with minimal heteroatoms
    # Assuming more than 12 carbon atoms considering the variances in chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, "Too few carbons to be a fatty acid"

    # Find carboxylic acid group typically at the end of the chain
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Ensure the molecule is not overly complex beyond typical epoxy fatty acid features
    # Avoid presence of additional rings, heterocycles, or large functional groups that are atypical
    # for a simple fatty acid structure

    # Examine bond count or complexity for confounding structures
    ring_systems = Chem.GetSymmSSSR(mol)
    if len(ring_systems) > 1:
        return False, "Too many rings for typical epoxy fatty acid"

    # Ensure stereochemistry of epoxide if needed, assuming it enhances specificity
    # This is commented out as RDKit handles stereochemistry incompleteness

    return True, "Contains an epoxide ring and characteristics of a epoxy fatty acid (long carbon chain and carboxylic acid group)"