"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: CHEBI:27310 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    An N-acyl-L-alpha-amino acid is an L-alpha-amino acid carrying an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N-acyl-L-alpha-amino acid backbone
    n_acyl_l_alpha_amino_acid_pattern = Chem.MolFromSmarts("[C@@H](N[C@@H](C(=O)))(C(=O)[O;H,-])[C@@H]")
    if not mol.HasSubstructMatch(n_acyl_l_alpha_amino_acid_pattern):
        return False, "No N-acyl-L-alpha-amino acid backbone found"

    # Check for stereochemistry (L-configuration)
    try:
        mol = AllChem.AssignAtomChiralTagsFromStructure(mol)
    except Exception as e:
        return False, f"Error in assigning stereochemistry: {str(e)}"

    chiral_centers = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED]
    if len(chiral_centers) < 2:
        return False, "Not enough chiral centers to determine L-configuration"

    for idx in chiral_centers:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetHybridization() == Chem.HybridizationType.SP3 and atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return False, "Stereochemistry not consistent with L-configuration"

    return True, "Contains N-acyl-L-alpha-amino acid backbone"