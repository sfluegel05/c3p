"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: CHEBI:25710 alkanethiol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for -SH group
    sh_pattern = Chem.MolFromSmarts("[SH]")
    if not mol.HasSubstructMatch(sh_pattern):
        return False, "No sulfanyl (-SH) group found"
    
    # Look for linear or branched alkyl/alkenyl chains attached to sulfur
    alkyl_pattern = Chem.MolFromSmarts("[SH]~[CX4,CX3]~*")
    if not mol.HasSubstructMatch(alkyl_pattern):
        return False, "No alkyl/alkenyl chain attached to sulfur"

    # Exclude peptides and amino acids
    aa_pattern = Chem.MolFromSmarts("[NX3H2,NX4H3+]")
    if mol.HasSubstructMatch(aa_pattern):
        return False, "Molecule appears to be a peptide or amino acid"

    return True, "Contains sulfanyl (-SH) group attached to an alkyl/alkenyl chain"