"""
Classifies: CHEBI:26607 saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    Saturated fatty acids have a terminal carboxylic acid group (-COOH), only carbon-carbon single bonds, and a long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the terminal carboxylic acid group (CX3(=O)O)
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) != 1:
        return False, f"Found {len(acid_matches)} carboxylic acid groups, must be exactly 1"

    # Get indices of the atoms belonging to the acid group
    acid_atom_indices = mol.GetSubstructMatch(acid_pattern)
    
    # Check for absence of C=C or C#C bonds
    unsaturated_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    if mol.HasSubstructMatch(unsaturated_pattern):
          return False, "Contains carbon-carbon double bonds"
    triple_pattern = Chem.MolFromSmarts("[CX2]#[CX2]")
    if mol.HasSubstructMatch(triple_pattern):
          return False, "Contains carbon-carbon triple bonds"
    
    # Check for minimum carbon chain length (arbitrarily set to 4, might need to be increased in future)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Too few carbons to be a fatty acid"
    
    # Check for only C, H, O, or 2H atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1 and atom.GetAtomicNum() != 8 and atom.GetAtomicNum() != 1:
          if atom.GetIsotope() != 2:
            return False, "Contains heteroatoms other than C, H, O or D"
        
        #Check that all oxygens belong to the carboxylic acid
        if atom.GetAtomicNum() == 8: #It's an oxygen atom
           if atom.GetIdx() not in acid_atom_indices:
               return False, "Contains O atoms outside of the terminal carboxyl group"
    
    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
          return False, "Contains a ring, fatty acids are acyclic"

    #Ensure that the carboxyl group is terminal
    
    carbonyl_carbon_index = acid_atom_indices[0]
    carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_index)

    neighbors = [atom.GetIdx() for atom in carbonyl_carbon.GetNeighbors()]
    if len(neighbors) > 2:
        return False, "The carbonyl group is not terminal"
    
    return True, "Contains a terminal carboxylic acid group, only carbon-carbon single bonds, no rings and a long carbon chain"