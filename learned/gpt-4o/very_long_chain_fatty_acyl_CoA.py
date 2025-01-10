"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    A very long-chain fatty acyl-CoA has a fatty acyl group with a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for coenzyme A substructure (simplified as certain patterns involving phosphates)
    coa_pattern = Chem.MolFromSmarts("COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](O)[C@H](O)[C@@H]1OP(O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A backbone found"

    # Find the longest carbon chain
    carbon_chain = max((len(list(filter(lambda x: x.GetAtomicNum() == 6, Chem.GetMolFrags(mol, asMols=True)[i].GetAtoms()))) for i in range(len(Chem.GetMolFrags(mol, asMols=True)))), default=0)

    if carbon_chain <= 22:
        return False, f"Carbon chain length is {carbon_chain}, not greater than C22"

    return True, f"Contains coenzyme A and a fatty acyl chain length of {carbon_chain} which is greater than C22"