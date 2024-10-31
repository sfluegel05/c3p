from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_D_tryptophan_derivative(smiles: str):
    """
    Determines if a molecule is a D-tryptophan derivative.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a D-tryptophan derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for indole core structure
    indole_pattern = Chem.MolFromSmarts('c1ccc2[nH]ccc2c1')
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole core structure found"

    # Check for amino acid backbone with correct stereochemistry
    # [C@H] ensures D-stereochemistry (clockwise)
    d_aa_pattern = Chem.MolFromSmarts('[C@H](CN)C(=O)[OH,O-,N]')
    if not mol.HasSubstructMatch(d_aa_pattern):
        return False, "No D-amino acid backbone found"

    # Check connection between indole and amino acid
    tryptophan_pattern = Chem.MolFromSmarts('c1ccc2[nH]ccc2c1CC[C@H](N)C(=O)[OH,O-,N]')
    if not mol.HasSubstructMatch(tryptophan_pattern):
        return False, "Indole not properly connected to amino acid backbone"

    # Get modification details
    modifications = []
    
    # Check for N-terminal modifications
    n_acetyl = Chem.MolFromSmarts('CC(=O)N[C@H]')
    if mol.HasSubstructMatch(n_acetyl):
        modifications.append("N-acetylated")
    
    n_malonyl = Chem.MolFromSmarts('OC(=O)CC(=O)N[C@H]')
    if mol.HasSubstructMatch(n_malonyl):
        modifications.append("N-malonylated")

    # Check for C-terminal modifications
    if "C(=O)N" in smiles:
        modifications.append("C-terminal amide")
        
    # Check for ring substitutions
    chloro_pattern = Chem.MolFromSmarts('Clc1cc2[nH]ccc2cc1')
    if mol.HasSubstructMatch(chloro_pattern):
        modifications.append("chlorinated")
        
    hydroxy_pattern = Chem.MolFromSmarts('Oc1cc2[nH]ccc2cc1')
    if mol.HasSubstructMatch(hydroxy_pattern):
        modifications.append("hydroxylated")

    if not modifications:
        modifications.append("unmodified D-tryptophan core")
        
    return True, "D-tryptophan derivative with modifications: " + ", ".join(modifications)
# Pr=None
# Recall=None