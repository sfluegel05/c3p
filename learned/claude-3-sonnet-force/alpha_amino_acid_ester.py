"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: CHEBI:38791 alpha-amino acid ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is an amino acid ester derivative obtained by the formal condensation of an alpha-amino acid with an alcohol.

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
    
    # Find alpha-amino acid backbone
    amino_acid_pattern = Chem.MolFromSmarts("[C](=[O])([N])([C])")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    # Check for alpha-amino acid backbone
    if not amino_acid_matches:
        return False, "No alpha-amino acid backbone found"
    
    # Find ester groups
    ester_pattern = Chem.MolFromSmarts("[C]([C](=[O])([O]))")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check for ester group
    if not ester_matches:
        return False, "No ester group found"
    
    # Check if ester is attached to alpha carbon of amino acid
    for amino_acid_match in amino_acid_matches:
        amino_acid_alpha_carbon = amino_acid_match[2]
        for ester_match in ester_matches:
            ester_carbon = ester_match[0]
            if mol.GetBondBetweenAtoms(amino_acid_alpha_carbon, ester_carbon).GetBondType() == BondType.SINGLE:
                return True, "Contains an alpha-amino acid backbone with an ester group attached to the alpha carbon"
    
    return False, "Ester group not attached to alpha carbon of amino acid backbone"