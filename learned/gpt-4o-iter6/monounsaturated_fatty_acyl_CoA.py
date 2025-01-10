"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    A monounsaturated fatty acyl-CoA contains a fatty acyl chain with one carbon-carbon double bond 
    and a CoA group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the CoA group substructure - includes larger moiety components
    coa_pattern = Chem.MolFromSmarts("COP(O)(=O)OP(O)(=O)OC1[C@@H](O)C(C)(C)N2C=NC3=C(N)N=CN=C3N2[C@H]1O[C@@H]2[C@H](O)[C@@H](O)[C@H](O2)OP(O)(O)=O") 
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA group found"

    # Define a pattern for carbon-carbon double bonds in aliphatic chains
    double_bond_pattern = Chem.MolFromSmarts("[CH2]=[CH]")  # Aliphatic C=C bond
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)

    # Research chain connections beyond just esters
    fatty_acid_chain_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)C")
    if not mol.HasSubstructMatch(fatty_acid_chain_pattern):
        return False, "No valid fatty acyl chain found"
    
    # Ensure there's exactly one carbon-carbon double bond
    num_double_bonds = len(double_bond_matches)
    if num_double_bonds != 1:
        return False, f"Found {num_double_bonds} double bonds in aliphatic chain, need exactly 1"

    return True, "Contains CoA and a fatty acyl chain with one carbon-carbon double bond"