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

    # Check for N-acyl-L-alpha-amino acid pattern.
    # This SMARTS pattern looks for an L-alpha-amino acid with an N-acyl group.
    # It looks for C(=O)-[C@H]-N-C(=O) where the C(=O) is attached to the N of the amino acid.
    n_acyl_amino_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NH1X2][C@H](C(=O)O)")
    if not mol.HasSubstructMatch(n_acyl_amino_acid_pattern):
         return False, "Molecule does not match N-acyl-L-alpha-amino acid pattern"

    return True, "Molecule is an N-acyl-L-alpha-amino acid"