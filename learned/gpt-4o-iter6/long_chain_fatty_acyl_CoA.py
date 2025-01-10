"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA includes a long-chain fatty acid esterified with the thiol group of coenzyme A.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the Coenzyme A motif, given its typical fragment features
    coa_pattern = Chem.MolFromSmarts("OP(O)(=O)O[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(=O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    # Define thioester linkage pattern (thioester should be near Coenzyme A)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Count carbon atoms in the chain (located before the ester group)
    for match in thioester_matches:
        carbon_chain_atoms = []  # Atoms constituting the fatty chain
        for bond in mol.GetAtomWithIdx(match[0]).GetNeighbors():
            # Trace unfettered linear carbon chains
            if bond.GetAtomicNum() == 6 and bond.GetSymbol() == 'C':
                carbon_chain_atoms.append(bond.GetIdx())
        
        carbon_count = len(carbon_chain_atoms)
        
        if 13 <= carbon_count <= 22:
            return True, f"Valid long-chain fatty acyl-CoA with {carbon_count} C atoms in chain"

    return False, "Carbon chain length not in C13-C22 for long-chain fatty acid"