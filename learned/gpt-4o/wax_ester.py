"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is a fatty acid ester resulting from the condensation of a fatty acid with a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ester functional group pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Check if the molecule contains an ester group
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"

    # Use SMARTS to check for long carbon chains, at least 10 carbons long
    long_chain_pattern = Chem.MolFromSmarts("[C][C][C][C][C][C][C][C][C][C]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long carbon chain (10+ carbons) found"
    
    # Estimate the molecular weight - wax esters are typically large
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for wax ester"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # A rough estimate: wax esters usually have many carbon atoms (e.g., more than 20)
    if c_count < 20:
        return False, "Too few carbon atoms for a wax ester"

    # Check for possible unsaturations in the alkyl chains
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return True, "Contains ester group with long carbon chains - likely a saturated wax ester"
    else:
        return True, "Contains ester group with long carbon chains including double bonds - likely an unsaturated wax ester"