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

    # Define a flexible ester pattern
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    
    # Check if the ester pattern is present in the molecule
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester groups found"

    # Define a 14-carbon chain pattern (not necessarily linear)
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]")
    
    # Check if the 14-carbon chain pattern is present in the molecule
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if len(carbon_chain_matches) < 1:
        return False, "No 14-carbon chain found"

    # Check molecular weight to ensure it's reasonable for a tetradecanoate ester
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:  # Adjusted to be more lenient
        return False, "Molecular weight too low for tetradecanoate ester"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 14:
        return False, "Too few carbons for tetradecanoate ester"
    if o_count < 2:
        return False, "Must have at least 2 oxygens (ester group)"

    return True, "Contains a tetradecanoate ester group (14-carbon chain with ester linkage)"