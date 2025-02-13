"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: CHEBI:64514 2-enoyl-CoA
An unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CoA scaffold (flexible pattern)
    coa_pattern = Chem.MolFromSmarts("[N&R]1C=NC2=NC=NC(=N1)N2OC[C@H]3O[C@@H]([C@H](O)[C@@H]3OP(O)(O)=O)COP(O)(=O)OCC(=O)NCCC(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA scaffold not found"

    # Look for double bond in S-acyl chain
    acyl_chain = Chem.MolFromSmarts("C=C")
    acyl_matches = mol.GetSubstructMatches(acyl_chain)

    # Check if double bond is part of the S-acyl group
    for match in acyl_matches:
        acyl_atoms = set()
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in match or bond.GetEndAtomIdx() in match:
                acyl_atoms.add(bond.GetBeginAtomIdx())
                acyl_atoms.add(bond.GetEndAtomIdx())

        acyl_atoms_coa = acyl_atoms.intersection(set(coa_pattern.GetSubstructMatch(mol)))
        if acyl_atoms_coa:
            break
    else:
        return False, "Double bond not part of S-acyl group"

    # Check for double bond stereochemistry (optional)
    cis_trans_bond = mol.GetBondWithIdx(mol.GetBonds()[list(acyl_atoms)[0].GetBondIdx(list(acyl_atoms)[1])])
    if cis_trans_bond.GetStereo() == Chem.BondStereo.STEREOE:
        stereochemistry = "trans"
    elif cis_trans_bond.GetStereo() == Chem.BondStereo.STEREOZ:
        stereochemistry = "cis"
    else:
        stereochemistry = "unknown"

    return True, f"Contains CoA scaffold with a {stereochemistry} double bond in the S-acyl group"