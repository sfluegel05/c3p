"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester has an alpha-amino acid core where the carboxyl group
    is replaced by an ester group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for alpha carbon with amino and carboxyl/ester groups (more specific)
    # This pattern explicitly requires the amino group and the ester group are attached to the same carbon
    # and the ester oxygen to a carbonyl carbon
    alpha_amino_acid_ester_pattern = Chem.MolFromSmarts("[CX4H0-1]([NX3H0-2])([CX3](=[OX1])[OX2][#6])")
    
    #check that the pattern is matched
    matches = mol.GetSubstructMatches(alpha_amino_acid_ester_pattern)

    if not matches:
      return False, "Molecule does not contain an alpha-amino acid with an ester group"
    
    for match in matches:
        #check if the carbon of the carboxylate has any other ester bonds attached
        carbonyl_carbon_idx = match[2]
        carbonyl_carbon_atom = mol.GetAtomWithIdx(carbonyl_carbon_idx)
        
        ester_bonds_count = 0
        for bond in carbonyl_carbon_atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(carbonyl_carbon_atom).GetAtomicNum()==8:
                ester_bonds_count +=1
        if ester_bonds_count != 1:
            return False, "Carbonyl carbon should not be part of another ester"
    return True, "Molecule contains an alpha-amino acid ester."