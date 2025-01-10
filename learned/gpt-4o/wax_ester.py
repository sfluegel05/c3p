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

    # Ensure no additional ester groups indicating complex structures (e.g., triglycerides)
    if mol.GetSubstructMatches(ester_pattern) and len(mol.GetSubstructMatches(ester_pattern)) > 1:
        return False, "Multiple ester groups, structure too complex"

    # Use SMARTS to check for potential fatty acid-alcohol ester with at least two separate long hydrocarbon chains
    long_chain_pattern = Chem.MolFromSmarts("[CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0]")
    chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    
    if len(chain_matches) < 2:
        return False, "Insufficient long carbon chains"

    # Check if structure is simple enough, avoiding branching and rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings, unlikely to be a simple wax ester"

    if any(atom.GetDegree() > 3 for atom in mol.GetAtoms()):
        return False, "Complex branching, unlikely to be wax ester"

    # Check the length of carbon chains on either side of the ester group
    ester_bond = mol.GetSubstructMatch(ester_pattern)
    if not ester_bond:
        return False, "Ester group not matching expected structure"
    
    c_below_chains = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1 and atom.GetIdx() < ester_bond[0]])
    c_above_chains = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1 and atom.GetIdx() > ester_bond[-1]])
    
    if c_below_chains < 10 or c_above_chains < 10:
        return False, "One or both chains are too short"

    # Verify the molecular weight - wax esters are typically large
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 450:
        return False, "Molecular weight below typical for wax ester"

    # Check for unsaturation
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    has_double_bond = mol.HasSubstructMatch(double_bond_pattern)
    
    if has_double_bond:
        return True, "Contains single ester group with long carbon chains including double bonds - likely an unsaturated wax ester"
    else:
        return True, "Contains single ester group with long carbon chains - likely a saturated wax ester"