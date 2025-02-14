"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA is a coenzyme A that results from the formal condensation of the
    thiol group of coenzyme A with the carboxy group of any branched-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CoA backbone pattern
    coa_pattern = Chem.MolFromSmarts("C(C(C(=O)NCCC(=O)NCCS)O)(C)(CO[P@@](=O)(O)O[P@@](=O)(O)OC[C@H]1[C@H]([C@@H]([C@H](O1)N1C=NC2=C1N=CN=C2N)O)O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"

    # Look for branched fatty acid chain pattern
    branched_pattern = Chem.MolFromSmarts("[CX4H3]([CX4H3])([CX4H3])[CX3]")
    branched_matches = mol.GetSubstructMatches(branched_pattern)
    if not branched_matches:
        return False, "No branched fatty acid chain found"

    # Look for ester linkage between CoA and fatty acid
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester linkage found"

    # Check for connectivity between branched fatty acid and CoA via ester
    for branched_match in branched_matches:
        for ester_match in ester_matches:
            if any(atom_idx in branched_match for atom_idx in ester_match):
                return True, "Contains a branched fatty acid chain linked to CoA via an ester bond"

    return False, "Branched fatty acid not linked to CoA via ester bond"