"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:36639 branched-chain fatty acyl-CoA
A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A 
with the carboxy group of any branched-chain fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.

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

    # Look for CoA backbone pattern (more flexible pattern)
    coa_pattern = Chem.MolFromSmarts("[C@@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        coa_pattern = Chem.MolFromSmarts("[CH](OP(O)(=O)OP(O)(=O)OC[CH]1O[CH]([CH](O)[CH]1OP(O)(O)=O)n1cnc2c(N)ncnc12)")
        coa_matches = mol.GetSubstructMatches(coa_pattern)
        if not coa_matches:
            return False, "No CoA backbone found"

    # Look for branched fatty acid chains
    branched_chain_pattern = Chem.MolFromSmarts("[CX4H2][CX4H2]([CX4])[CX4H2]")
    branched_chain_matches = mol.GetSubstructMatches(branched_chain_pattern)
    if not branched_chain_matches:
        return False, "No branched fatty acid chain found"

    # Check for ester bond between CoA and fatty acid chain
    ester_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    ester_bond_matches = mol.GetSubstructMatches(ester_bond_pattern)
    if not ester_bond_matches:
        return False, "No ester bond found between CoA and fatty acid chain"

    # Check for long fatty acid chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Fatty acid chain too short"

    # Check for branching in the fatty acid chain
    for match in branched_chain_matches:
        atom1 = mol.GetAtomWithIdx(match[0])
        atom2 = mol.GetAtomWithIdx(match[1])
        atom3 = mol.GetAtomWithIdx(match[2])
        atom4 = mol.GetAtomWithIdx(match[3])

        # Check for branching at atom2 or atom3
        branched_at_2 = sum(1 for bond in atom2.GetBonds() if bond.GetBondType() == Chem.BondType.SINGLE) > 2
        branched_at_3 = sum(1 for bond in atom3.GetBonds() if bond.GetBondType() == Chem.BondType.SINGLE) > 2

        if branched_at_2 or branched_at_3:
            # Check for potential substitutions or functional groups
            substitutions = ['O', '=O', '=C']
            for atom in mol.GetAtoms():
                if atom.GetSymbol() in substitutions:
                    return True, "Contains CoA backbone with branched fatty acid chain and potential substitutions"

    return False, "No branching or relevant substitutions found in the fatty acid chain"