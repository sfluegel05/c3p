"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA
"""

from rdkit import Chem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA is a fatty acyl-CoA resulting from the condensation 
    of the thiol group of coenzyme A with the carboxy group of any branched-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    from rdkit import Chem

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find the thioester bond (C(=O)-S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Assume the first match is the thioester linkage
    thioester_match = thioester_matches[0]
    carbonyl_c_idx = thioester_match[0]
    sulfur_idx = thioester_match[2]
    bond = mol.GetBondBetweenAtoms(carbonyl_c_idx, sulfur_idx)
    if bond is None:
        return False, "No thioester bond found"

    # Break the molecule at the thioester bond to isolate the acyl chain
    bonds_to_break = [bond.GetIdx()]
    frags = Chem.FragmentOnBonds(mol, bonds_to_break, addDummies=False)
    # Get the fragments
    mol_frags = Chem.GetMolFrags(frags, asMols=True, sanitizeFrags=True)

    # Identify the acyl chain fragment (the one containing the carbonyl carbon)
    acyl_chain_frag = None
    for frag in mol_frags:
        frag_atom_indices = frag.GetSubstructMatch(frag)
        # Check if the fragment contains the carbonyl carbon
        if any(atom.GetAtomicNum() == 6 and atom.GetIdx() == carbonyl_c_idx for atom in frag.GetAtoms()):
            acyl_chain_frag = frag
            break
    if acyl_chain_frag is None:
        return False, "Could not isolate acyl chain fragment"

    # Check for rings in the acyl chain
    if acyl_chain_frag.GetRingInfo().NumRings() > 0:
        return False, "Acyl chain contains ring structures"

    # Check for heteroatoms other than oxygen (allowing for keto and hydroxyl groups)
    for atom in acyl_chain_frag.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in [1, 6, 8]:  # Allow only H, C, O
            return False, f"Acyl chain contains heteroatom {atom.GetSymbol()}"

    # Check for branching in the acyl chain
    branching_points = 0
    for atom in acyl_chain_frag.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            neighbor_carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if len(neighbor_carbons) > 2:
                branching_points += 1

    if branching_points == 0:
        return False, "Acyl chain is not branched"
    else:
        return True, f"Acyl chain is branched with {branching_points} branching point(s)"