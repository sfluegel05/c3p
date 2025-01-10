"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is an unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.

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

    # Identify thioester linkage: C(=O)-S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Assume first match is the thioester linkage
    carbonyl_c_idx = thioester_matches[0][0]
    sulfur_idx = thioester_matches[0][2]
    sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)

    # Check for CoA moiety connected to sulfur atom
    # Use a SMARTS pattern for the CoA fragment (pantetheine with diphosphate adenine)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(C(O)C1O)N2C=NC3=C2N=CN=C3N")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found connected to sulfur atom"

    # Identify the carbonyl carbon atom
    carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
    # Find the alpha carbon (next carbon in acyl chain)
    alpha_c = None
    for neighbor in carbonyl_c.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != sulfur_idx:
            alpha_c = neighbor
            break
    if alpha_c is None:
        return False, "Alpha carbon in acyl chain not found"

    # Check if alpha carbon is sp2 hybridized (double-bonded to beta carbon)
    if alpha_c.GetHybridization() != Chem.rdchem.HybridizationType.SP2:
        return False, "Alpha carbon is not sp2 hybridized"

    # Identify beta carbon (double-bonded to alpha carbon)
    beta_c = None
    for neighbor in alpha_c.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != carbonyl_c_idx:
            bond = mol.GetBondBetweenAtoms(alpha_c.GetIdx(), neighbor.GetIdx())
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                beta_c = neighbor
                break
    if beta_c is None:
        return False, "Beta carbon double-bonded to alpha carbon not found"

    # Ensure the double bond is not part of a ring
    if alpha_c.IsInRing() or beta_c.IsInRing():
        return False, "Double bond is part of a ring"

    # Check that the acyl chain continues beyond beta carbon
    chain_continues = False
    for neighbor in beta_c.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != alpha_c.GetIdx():
            chain_continues = True
            break
    if not chain_continues:
        return False, "Acyl chain does not continue beyond beta carbon"

    return True, "2-enoyl-CoA identified with double bond between positions 2 and 3 in acyl chain"