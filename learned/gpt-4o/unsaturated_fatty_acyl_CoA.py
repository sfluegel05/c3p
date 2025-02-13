"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
from rdkit import Chem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA results from the formal condensation of the thiol group of coenzyme A 
    with the carboxy group of an unsaturated fatty acid.
    
    Args:
    - smiles (str): SMILES string of the molecule
    
    Returns:
    - bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
    - str: Reason for classification
    """
    # Parse the SMILES string into an RDKit mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for thioester linkage between fatty acid chain and CoA: C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Improved pattern for Coenzyme A recognition, accounting for charged states
    coenzymeA_pattern = Chem.MolFromSmarts("NC(=O)[C@H](O)[C](C)(C)COP(=O)(O)O[C@H]")
    if not mol.HasSubstructMatch(coenzymeA_pattern):
        return False, "Coenzyme A moiety pattern not found"

    # Recognize unsaturated fatty acid chains with substantial chain length
    # Ensure there is an unsaturation (C=C bond) in the context of a long carbon chain
    unsaturated_fatty_acid_pattern = Chem.MolFromSmarts("C=C[CH2,CH]")
    if not mol.HasSubstructMatch(unsaturated_fatty_acid_pattern):
        return False, "No unsaturation (C=C bond) found in the fatty acid chain"

    # Check chain length, this ensures we're getting significant carbon chains
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCC")  # Minimum length example
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Fatty acid chain is too short"

    return True, "Molecule is an unsaturated fatty acyl-CoA, identified by the presence of CoA moiety with long unsaturated fatty acid chain"