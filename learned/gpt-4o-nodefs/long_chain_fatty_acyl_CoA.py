"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the Coenzyme A (CoA) structure pattern including adenine, ribose, and phosphates
    coa_smarts = "NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)([O-])=OOP(O)(=O)OC[C@@H]1O[C@H](CN2C=NC3=C2N=CN=C3N)C(O)C1O"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing or incomplete Coenzyme A structure"
    
    # Define the thiol ester linkage
    thiol_ester_smarts = "C(=O)SCCNC(=O)"
    thiol_ester_pattern = Chem.MolFromSmarts(thiol_ester_smarts)
    if not mol.HasSubstructMatch(thiol_ester_pattern):
        return False, "Missing proper thiol ester linkage"

    # Recognize a long hydrocarbon chain
    long_chain_exists = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            chain_length = 1
            neighbors = atom.GetNeighbors()
            while len(neighbors) == 1 and neighbors[0].GetAtomicNum() == 6:
                chain_length += 1
                neighbors = neighbors[0].GetNeighbors()
            if chain_length >= 16:
                long_chain_exists = True
                break
    if not long_chain_exists:
        return False, "No sufficiently long hydrocarbon chain found"

    # Check if there are acceptable numbers of double bonds
    db_pattern = Chem.MolFromSmarts("C=C")
    db_matches = mol.GetSubstructMatches(db_pattern)
    if len(db_matches) > 6:
        return False, "Too many double bonds for typical long-chain fatty acyl-CoA"
    
    return True, "Matches long-chain fatty acyl-CoA features"