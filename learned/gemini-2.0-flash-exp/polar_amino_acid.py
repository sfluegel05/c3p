"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid has a side chain capable of forming hydrogen bonds.

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

    # Check for the amino acid backbone pattern. Relaxed nitrogen hydrogens, allowing for zwitterions and protonations
    backbone_pattern = Chem.MolFromSmarts("[NX3]-[CX4](-[CX3](=[OX1])-[OX2;H1])")
    if not mol.HasSubstructMatch(backbone_pattern):
         return False, "Molecule does not have the amino acid backbone"

    # Find the alpha carbon (the one between the N and carbonyl C)
    alpha_carbon_pattern = Chem.MolFromSmarts("[NX3]-[CX4](-[CX3](=[OX1])-[OX2;H1])")
    matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    if not matches:
      return False, "Could not find alpha carbon"
    alpha_carbon_idx = matches[0][1] #index 1 corresponds to the alpha carbon
   
    
    # Check for polar side chain patterns attached to sidechains
    polar_side_chain_patterns = [
        Chem.MolFromSmarts("[CX4]-[OX2;H1]"),  # hydroxyl (-OH)
        Chem.MolFromSmarts("[CX4]-[NX3;H2,H1]"),  # amino (-NH2)
        Chem.MolFromSmarts("[CX4]-C(=[OX1])N"),  # amide (-CONH2)
        Chem.MolFromSmarts("[CX4]-C(=[OX1])[OX2;H1]"),  # carboxylic acid (-COOH)
        Chem.MolFromSmarts("[CX4]-[SX2;H1]"), #thiol (-SH)
        Chem.MolFromSmarts("[CX4]-c[nH]c"),   #imidazole
        Chem.MolFromSmarts("[CX4]-c1ccccc1[OH]"), #phenol
        Chem.MolFromSmarts("[CX4]-c1ccncc1N") #pyridine
    ]
    
    has_polar_group = False
    for pattern in polar_side_chain_patterns:
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            polar_carbon_idx = match[0]
            if mol.GetAtomWithIdx(polar_carbon_idx).GetIdx() != alpha_carbon_idx:
              has_polar_group = True
              break
        if has_polar_group:
            break
           
    if not has_polar_group:
        return False, "No polar side chain capable of forming hydrogen bonds found."

    #check that there are no peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3]-[CX3](=[OX1])-[NX3]")
    if mol.HasSubstructMatch(peptide_bond_pattern):
       return False, "Molecule contains a peptide bond and therefore is not a single amino acid."
    

    return True, "Molecule has amino acid backbone with polar side chain"