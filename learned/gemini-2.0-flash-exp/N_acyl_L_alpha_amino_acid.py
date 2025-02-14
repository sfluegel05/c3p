"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    #Check for stereochemistry
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    for atom_index, chirality_tag in chiral_centers:
       if chirality_tag == '?':
        continue
       atom = mol.GetAtomWithIdx(atom_index)
       if atom.GetSymbol() == 'C' and atom.HasProp('_Chirality'):
           if atom.GetProp('_Chirality') == 'S': # All chiral carbons are L
            return False, f"Molecule does not have L configuration on all chiral carbons. Atom {atom_index} is S"
           elif not atom.GetProp('_Chirality') == 'R':
            return False, f"Molecule does not have L configuration on all chiral carbons. Atom {atom_index} is not R"
           

    # Check for L-alpha-amino acid core pattern
    amino_acid_pattern = Chem.MolFromSmarts("[NX2][C@H](C(=O)O)")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Molecule does not contain L-alpha-amino acid core."
    
    # Check for N-acyl group (-C(=O)-) attached to any N atom.
    acyl_group_pattern = Chem.MolFromSmarts("[NX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(acyl_group_pattern):
         return False, "Molecule does not contain an N-acyl group"

    return True, "Molecule is an N-acyl-L-alpha-amino acid"