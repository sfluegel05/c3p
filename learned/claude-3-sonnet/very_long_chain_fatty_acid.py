"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:36216 very long-chain fatty acid

A fatty acid which has a chain length greater than C22. Very long-chain fatty acids
which have a chain length greater than C27 are also known as ultra-long-chain fatty acids.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    has_acid = any(atom.GetSymbol() == 'O' and sum(mol.GetAtomWithIdx(i).GetTotalNumHs() for i in atom.GetNeighbors()) == 1
                   for atom in mol.GetAtoms())
    if not has_acid:
        return False, "Does not contain a carboxylic acid group"
    
    # Get the longest carbon chain
    chain = Chem.Lipinski.GetLongestChain(mol)
    if chain is None:
        return False, "No carbon chain found"

    # Check chain length
    chain_length = len(chain)
    if chain_length <= 22:
        return False, f"Carbon chain length ({chain_length}) too short"
    
    # Check for unsaturations
    unsaturations = Chem.Lipinski.GetNumUnsaturatedRings(mol) + Chem.Lipinski.GetNumAlkeneRings(mol)
    if unsaturations > 4:
        return False, "More than 4 unsaturations"
    
    # Check for allowed atoms
    allowed_atoms = {'C', 'H', 'O'}
    atoms = set(atom.GetSymbol() for atom in mol.GetAtoms())
    if not atoms.issubset(allowed_atoms):
        return False, "Contains disallowed atoms"
    
    # Check for branching
    if len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C]([C])([C])[C]'))) > 0:
        return False, "Contains branching"
    
    # Check for cycles
    if any(len(ring) > 8 for ring in mol.GetRingInfo().AtomRings()):
        return False, "Contains large cycles"
    
    # Check for substituents
    if any(atom.GetTotalNumHs(False) > 1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O'):
        return False, "Contains hydroxyl or ether substituents"
    if any(atom.GetTotalNumHs(False) == 0 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O'):
        return False, "Contains carbonyl or ester substituents"
    
    # Check stereochemistry
    if mol.GetBondBetweenAtoms(0, 1).GetStereo() != Chem.BondStereo.STEREONONE:
        return False, "Contains stereochemistry"
    
    return True, "Meets all criteria for a very long-chain fatty acid"