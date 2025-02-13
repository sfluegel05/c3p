"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: CHEBI:44243 alpha-amino acid ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is the amino acid ester derivative obtained by the formal
    condensation of an alpha-amino acid with an alcohol.

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
    
    # Look for alpha-amino acid backbone
    aa_backbone_pattern = Chem.MolFromSmarts("[NX3H2,NX4H3+0][CX4H]([CH0-1])(C(=O)[OX1H0-,OX2H1])[CH0-1]")
    if not mol.HasSubstructMatch(aa_backbone_pattern):
        return False, "No alpha-amino acid backbone found"
    
    # Look for ester group
    ester_pattern = Chem.MolFromSmarts("C(=O)[OX2H0]C")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Check if ester is connected to alpha carbon
    aa_ester_pattern = Chem.MolFromSmarts("[NX3H2,NX4H3+0][CX4H]([CH0-1])(C(=O)[OX1H0-,OX2H1])[CH0-1]C(=O)[OX2H0]C")
    if not mol.HasSubstructMatch(aa_ester_pattern):
        return False, "Ester group not connected to alpha-amino acid backbone"
    
    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 500:
        return False, "Molecular weight outside typical range for alpha-amino acid esters"
    
    n_atoms = mol.GetNumAtoms()
    if n_atoms < 10 or n_atoms > 60:
        return False, "Number of atoms outside typical range for alpha-amino acid esters"
    
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2 or n_rotatable > 20:
        return False, "Number of rotatable bonds outside typical range for alpha-amino acid esters"
    
    return True, "Molecule contains an alpha-amino acid backbone with an ester group attached"