"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: CHEBI:27279 vitamin D
Any member of a group of fat-soluble hydroxy seco-steroids that exhibit biological activity against vitamin D deficiency. 
Vitamin D can be obtained from sun exposure, food and supplements and is biologically inactive and converted into 
the biologically active calcitriol via double hydroxylation in the body.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for seco-steroid backbone
    seco_steroid_pattern = Chem.MolFromSmarts("[C@]12[C@H](CCC[C@@H]1[C@H](C)CCCC(C)(C)O)[C@@H](CC[C@@]3([C@H]2CCC4=CC(=O)CC[C@@]34C)C)O"
    if not mol.HasSubstructMatch(seco_steroid_pattern):
        return False, "No seco-steroid backbone found"
    
    # Check for cis-triene system
    cis_triene_pattern = Chem.MolFromSmarts("/C=C/C=C/C=C/"
    if not mol.HasSubstructMatch(cis_triene_pattern):
        return False, "No cis-triene system found"
    
    # Check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Less than 2 hydroxyl groups found"
    
    # Check molecular weight range (typically 300-500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 500:
        return False, "Molecular weight outside typical range for vitamin D"
    
    # Check for additional structural features
    additional_features_pattern = Chem.MolFromSmarts("[C@@]1(CCC[C@@H](C1)O)[C@H](C)CCCC(C)(C)O"
    if not mol.HasSubstructMatch(additional_features_pattern):
        return False, "Missing additional structural features of vitamin D"
    
    return True, "Molecule exhibits structural features of vitamin D"