"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: CHEBI:36601 tetradecanoate ester
A fatty acid ester obtained by condensation of the carboxy group of tetradecanoic acid 
(also known as myristic acid) with a hydroxy group of an alcohol or phenol.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for tetradecanoic acid (myristic acid) substructure
    myristic_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O")
    myristic_acid_matches = mol.GetSubstructMatches(myristic_acid_pattern)
    if not myristic_acid_matches:
        return False, "No myristic acid substructure found"
    
    # Look for ester bonds
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester bonds found"
    
    # Check for alcohol/phenol group
    alcohol_pattern = Chem.MolFromSmarts("O[C;H2]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    phenol_pattern = Chem.MolFromSmarts("Oc")
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if not alcohol_matches and not phenol_matches:
        return False, "No alcohol or phenol groups found"
    
    # Check molecular weight - tetradecanoate esters typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for tetradecanoate ester"
    
    return True, "Contains myristic acid substructure connected via ester bond to an alcohol or phenol group"