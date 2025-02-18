"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is an ester of tetradecanoic acid (myristic acid) with an alcohol or phenol.

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

    # Check for the tetradecanoyl group attached to the *carbon* of an ester group.
    tetradecanoyl_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH3]")
    if not mol.HasSubstructMatch(tetradecanoyl_ester_pattern):
        return False, "No tetradecanoyl group attached to the carbon of an ester group found"

    # Verify minimum number of carbons (14 in tetradecanoyl and at least one more carbon in the alcohol chain)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
         return False, "Too few carbons to be a tetradecanoate ester"
    
    # Verify that molecular weight is within the appropriate range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for tetradecanoate ester"

    return True, "Contains tetradecanoyl group in an ester linkage"