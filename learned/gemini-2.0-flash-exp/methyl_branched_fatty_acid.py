"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid has a carboxylic acid group, a hydrocarbon chain,
    and only methyl branches.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

     # Check for multiple bonds outside carbonyl
    multiple_bond = Chem.MolFromSmarts("[!C](=[C,O,N,S])")
    if mol.HasSubstructMatch(multiple_bond):
        return False, "Molecule contains double or triple bonds"

    # Check for ring structures - not allowed
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[R]")):
        return False, "Molecule contains rings"

    #Check for other atoms except C,H,O
    other_atoms = Chem.MolFromSmarts("[!#1!#6!#8]")
    if mol.HasSubstructMatch(other_atoms):
         return False, "Molecule contains non C, H, or O atoms"

    # Ensure there are at least 4 carbons in the chain
    # excluding methyl carbons.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    methyl_group_pattern = Chem.MolFromSmarts("[CX4H3]")
    methyl_matches = mol.GetSubstructMatches(methyl_group_pattern)
    num_methyl_carbons = len(methyl_matches)
    num_carbons = len(carbon_atoms)
    if num_carbons - num_methyl_carbons < 4:
      return False, "Too few carbons in chain to be a fatty acid"
    
    # Check for minimum carbon chain length. (5 carbons, excluding those in branches)
    if num_carbons < 5:
        return False, "Too few carbons to be a fatty acid"

    #Find carbon connected to carboxyl group
    chain_carbon_pattern = Chem.MolFromSmarts("CC(=O)O")
    chain_carbon_matches = mol.GetSubstructMatches(chain_carbon_pattern)
    if not chain_carbon_matches:
        return False, "Carboxylic acid not attached to carbon chain"
    
    # Ensure all branches are methyl branches
    # Look for carbon atoms with more than 2 neighbors, and make sure they are connected only to CH3s, H or C.
    # A methyl branch is defined as a carbon atom with three non-hydrogen neighbours, one of which is a methyl (CH3).
    branch_pattern = Chem.MolFromSmarts("[CX4]([CH3])([#1,#6])([#1,#6])")
    methyl_branches = mol.GetSubstructMatches(branch_pattern)

   # Check for non-methyl branches. A non-methyl branch is a carbon attached to any non-hydrogen atom or carbon.
    non_methyl_branch_pattern = Chem.MolFromSmarts("[CX4]([!#1,#6])([!#1,#6])([!#1,#6])")
    if mol.HasSubstructMatch(non_methyl_branch_pattern):
        return False, "Molecule contains non-methyl branches"

    return True, "Methyl-branched fatty acid"