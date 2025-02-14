"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid has a side chain capable of forming one or more hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): True if molecule is a polar amino acid, False otherwise
                         and the reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the amino acid backbone pattern
    # Relaxed Nitrogen and Oxygen, also taking into account Deuterium.
    backbone_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]-[CX4](-[CX3](=[OX1])-[OX2;H1,H0])")
    if not mol.HasSubstructMatch(backbone_pattern):
         return False, "Molecule does not have the amino acid backbone"

    # Find the alpha carbon (the one between the N and carbonyl C)
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]-[CX4](-[CX3](=[OX1])-[OX2;H1,H0])")
    matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    if not matches:
      return False, "Could not find alpha carbon"
    alpha_carbon_idx = matches[0][1] #index 1 corresponds to the alpha carbon
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)

    #Get the neighbours of the alpha carbon
    neighbors = [atom for atom in alpha_carbon.GetNeighbors() if atom.GetIdx() != matches[0][0] and atom.GetIdx() != matches[0][2]]

    #check if any neighbour is polar
    has_polar_group = False
    for neighbor in neighbors:
      atomic_num = neighbor.GetAtomicNum()
      if atomic_num == 8:
        has_polar_group = True #O
        break
      if atomic_num == 7:
        has_polar_group = True #N
        break
      if atomic_num == 16:
          has_polar_group=True #S
          break
    if not has_polar_group:
        return False, "No polar side chain capable of forming hydrogen bonds found."
    
     #check that there are no peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3]-[CX3](=[OX1])-[NX3]")
    if mol.HasSubstructMatch(peptide_bond_pattern):
       return False, "Molecule contains a peptide bond and therefore is not a single amino acid."


    return True, "Molecule has amino acid backbone with polar side chain"