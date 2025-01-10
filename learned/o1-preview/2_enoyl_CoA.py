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

    # Check for CoA moiety
    # CoA contains an adenosine moiety with a ribose ring and adenine base
    adenosine_smarts = "n1cnc2c1ncnc2N[C@H]3O[C@H]([C@H](O)[C@@H]3O)CO"
    adenosine_pattern = Chem.MolFromSmarts(adenosine_smarts)
    if not mol.HasSubstructMatch(adenosine_pattern):
        return False, "Adenosine moiety of CoA not found"

    # Check for phosphate groups (simplified as phosphorus atoms with four oxygens)
    phosphate_smarts = "P(=O)(O)O"
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 3:
        return False, "Less than 3 phosphate groups found in CoA moiety"

    # Check for thioester linkage: C(=O)S
    thioester_smarts = "C(=O)S"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"

    # Check for double bond between carbons 2 and 3 in the acyl chain
    # Find the carbonyl carbon in thioester linkage
    carbonyl_smarts = "[#6](=O)S"
    carbonyl_pattern = Chem.MolFromSmarts(carbonyl_smarts)
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if not carbonyl_matches:
        return False, "Carbonyl carbon of thioester not found"

    # For each carbonyl carbon, check the connected chain for C=C at position 2
    for match in carbonyl_matches:
        carbonyl_c_idx = match[0]
        # Get the sulfur atom connected to carbonyl carbon
        sulfur_atom = mol.GetAtomWithIdx(match[1])
        # Get the carbon attached to sulfur (start of acyl chain)
        neighbors = [atom for atom in sulfur_atom.GetNeighbors() if atom.GetIdx() != carbonyl_c_idx]
        if not neighbors:
            continue
        acyl_chain_start = neighbors[0]
        # Check for double bond between carbons 2 and 3
        # Get neighbors of acyl_chain_start excluding sulfur
        acyl_neighbors = [atom for atom in acyl_chain_start.GetNeighbors() if atom.GetIdx() != sulfur_atom.GetIdx()]
        for atom in acyl_neighbors:
            # Check if there is a double bond to the next carbon
            bond = mol.GetBondBetweenAtoms(acyl_chain_start.GetIdx(), atom.GetIdx())
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and atom.GetAtomicNum() == 6:
                return True, "2-enoyl-CoA identified with double bond at position 2 in acyl chain"

    return False, "Double bond between positions 2 and 3 in acyl chain not found"