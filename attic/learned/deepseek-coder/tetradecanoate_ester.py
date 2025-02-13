"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is a fatty acid ester obtained by condensation of the carboxy group 
    of tetradecanoic acid (myristic acid) with a hydroxy group of an alcohol or phenol.

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

    # Define the tetradecanoate ester pattern: 14-carbon chain with ester linkage
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC(=O)O")
    
    # Check if the pattern is present in the molecule
    if not mol.HasSubstructMatch(tetradecanoate_pattern):
        return False, "No tetradecanoate ester group found"

    # Count the number of ester groups
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester groups found"

    # Check for the presence of a 14-carbon chain
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if len(carbon_chain_matches) < 1:
        return False, "No 14-carbon chain found"

    # Check molecular weight to ensure it's reasonable for a tetradecanoate ester
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:  # Tetradecanoate esters typically have MW > 200
        return False, "Molecular weight too low for tetradecanoate ester"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 14:
        return False, "Too few carbons for tetradecanoate ester"
    if o_count < 2:
        return False, "Must have at least 2 oxygens (ester group)"

    return True, "Contains a tetradecanoate ester group (14-carbon chain with ester linkage)"