"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the deprotonated form of a fatty acid, characterized by a carboxylate group (-COO-) at the end of a long hydrocarbon chain (>= 6 carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a terminal carboxylate group using SMARTS pattern
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "Molecule does not contain a terminal carboxylate group."
    
     #Find carboxylate matches
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) != 1:
        return False, "Molecule must contain exactly one carboxylate group"

    # Check that there are no other anions present (sulfonate/phosphate etc.)
    anion_pattern = Chem.MolFromSmarts("[!#6!#1][-,@]") # any non-carbon, non-hydrogen atom with negative charge
    anion_matches = mol.GetSubstructMatches(anion_pattern)
    if len(anion_matches) > 1 :
      return False, "Molecule contains more than one anionic group"

    # Check that there are no zwitterions present
    zwitterion_pattern = Chem.MolFromSmarts("[#7+]")
    zwitterion_matches = mol.GetSubstructMatches(zwitterion_pattern)
    if len(zwitterion_matches) > 0 :
      return False, "Molecule is a zwitterion"

    # Check for a carbon chain (>= 6 carbons) attached to carboxylate group
    # This pattern checks for a carbon connected to the carboxylate carbon, and then an additional 5 carbons
    chain_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[O-][CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    carbon_chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(carbon_chain_matches) == 0:
      return False, "Carboxylate group not attached to a carbon chain with length >= 6 carbons"

    # Check for rings and heteroatoms in the chain
    no_ring_hetero_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[O-][CX4,CX3](~[CX4,CX3])(~[CX4,CX3])(~[CX4,CX3])(~[CX4,CX3])(~[CX4,CX3])~[!#7,!#8,!#16,R]") # carbon chain with branching, no rings or heteroatoms
    no_ring_hetero_matches = mol.GetSubstructMatches(no_ring_hetero_pattern)
    if len(no_ring_hetero_matches) == 0:
       return False, "Carboxylate chain contains a ring or heteroatom"

    return True, "Contains a carboxylate group at the end of a carbon chain (>=6 carbons), consistent with a fatty acid anion."