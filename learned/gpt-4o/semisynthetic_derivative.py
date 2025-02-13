"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule could be classified as a semisynthetic derivative based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a semisynthetic derivative, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Assess complexity using molecular descriptors: MolWt and TPSA as proxies
    mol_wt = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)

    # Consider compounds above certain size/complexity as potential semisynthetic derivatives
    if mol_wt < 300 or tpsa < 40:
        return False, "Molecule seems too simple to be semisynthetic"
    
    # Look for typical semisynthetic modification groups
    acetyl_group = Chem.MolFromSmarts("CC(=O)O")
    if mol.HasSubstructMatch(acetyl_group):
        return True, "Contains acetyl group, indicates possible chemical modification"

    methylation = Chem.MolFromSmarts("C")
    methyl_count = len(mol.GetSubstructMatches(methylation))
    if methyl_count > 3:
        return True, "High degree of methylation, indicates possible semisynthetic modification"

    halogen_group = Chem.MolFromSmarts("[F,Cl,Br,I]")
    if mol.HasSubstructMatch(halogen_group):
        return True, "Contains halogen group, indicates possible chemical modification"
    
    amide_group = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_group):
        return True, "Contains amide group, common in semisynthetic derivatives"
    
    # If none of the inference patterns provide definitive classification
    return None, "Cannot reliably classify as semisynthetic derivative based on structure alone"