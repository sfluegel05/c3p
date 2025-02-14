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

    #Get the neighbours of the alpha carbon, excluding the backbone atoms
    neighbors = [atom for atom in alpha_carbon.GetNeighbors() if atom.GetIdx() != matches[0][0] and atom.GetIdx() != matches[0][2]]

    #Check for peptide bond
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3]-[CX3](=[OX1])-[NX3]")
    if mol.HasSubstructMatch(peptide_bond_pattern):
       return False, "Molecule contains a peptide bond and therefore is not a single amino acid."

    # Define SMARTS patterns for common polar groups in side chains
    polar_patterns = [
        Chem.MolFromSmarts("[OX2;H1]"), # hydroxyl
        Chem.MolFromSmarts("[NX3;H2]"), # primary amine
        Chem.MolFromSmarts("[NX3;H1][CX4]"), # secondary amine
        Chem.MolFromSmarts("[NX3;H0]([CX4])([CX4])"), # tertiary amine
        Chem.MolFromSmarts("C(=O)[NX2;H2]"), # amide
        Chem.MolFromSmarts("C(=O)[OX2;H1]"), # carboxylic acid
        Chem.MolFromSmarts("[SX2;H1]"),  #thiol
        Chem.MolFromSmarts("[NX2;H1]=[CX3][NX2;H2]"), #guanidine group
        Chem.MolFromSmarts("[NX3;H2][CX3]=[NX2]"), #guanidine group
        Chem.MolFromSmarts("[NX2;H1]c1[nX2;H1]cc[nX2]c1") #Imidazole group
        ]


    # Check if any of the polar patterns match atoms reachable from the alpha carbon, excluding the backbone atoms.
    has_polar_group = False
    for neighbor in neighbors:
        for pattern in polar_patterns:
            if mol.HasSubstructMatch(pattern, fromAtoms=[neighbor.GetIdx()]):
                has_polar_group = True
                break
        if has_polar_group:
            break

    if not has_polar_group:
            return False, "No polar side chain capable of forming hydrogen bonds found."

    return True, "Molecule has amino acid backbone with polar side chain"