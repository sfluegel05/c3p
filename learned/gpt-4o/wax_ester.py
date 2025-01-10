"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is a fatty acid ester resulting from the condensation of the carboxy group of a fatty acid with the alcoholic hydroxy group of a fatty alcohol.

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

    # Use SMARTS to check for potential fatty acid-alcohol ester with at least two long hydrocarbon chains
    long_chain_pattern1 = Chem.MolFromSmarts("[CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0]")  # Straight chain
    long_chain_pattern2 = Chem.MolFromSmarts("[CR0][CR0][CR0][CR0][CR0][CR0][CR0]=[CR0]")  # Chains can have some unsaturation
    chain_matches1 = mol.GetSubstructMatches(long_chain_pattern1)
    chain_matches2 = mol.GetSubstructMatches(long_chain_pattern2)

    if len(chain_matches1) + len(chain_matches2) < 2:
        return False, "Insufficient long carbon chains"

    # Verify that molecule is not complexly branched or cyclic
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings, less likely to be a wax ester"

    # Re-evaluate the molecular weight - wax esters are typically large
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight below typical for wax ester"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbon atoms for a wax ester"

    # Check for unsaturation
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    has_double_bond = mol.HasSubstructMatch(double_bond_pattern)
    
    if has_double_bond:
        return True, "Contains ester group with long carbon chains including double bonds - likely an unsaturated wax ester"
    else:
        return True, "Contains ester group with long carbon chains - likely a saturated wax ester"