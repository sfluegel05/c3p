"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA results from the condensation of the thiol group of coenzyme A
    with the carboxy group of any fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coenzyme A SMILES (from ChEBI: 15346)
    coa_smiles = "C[C@H](C(=O)NCCS)[C@H](O)COP(=O)(O)OP(=O)(O)OC[C@@H]1O[C@@H]" \
                 "([C@@H](O)[C@H]1O)n1cnc2c(N)ncnc12"
    
    # Create CoA molecule
    coa_mol = Chem.MolFromSmiles(coa_smiles)
    if coa_mol is None:
        return False, "Failed to create CoA molecule"

    # Modify CoA molecule to have attachment point at sulfur atom
    # Find the sulfur atom in CoA (should be only one sulfur atom)
    sulfur_idx = None
    for atom in coa_mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur atom
            sulfur_idx = atom.GetIdx()
            break
    if sulfur_idx is None:
        return False, "No sulfur atom found in CoA molecule"

    # Convert sulfur atom to a wildcard atom for matching thioester linkage
    sulfur_atom = coa_mol.GetAtomWithIdx(sulfur_idx)
    # Remove hydrogen atoms connected to sulfur
    sulfur_atom.SetNumExplicitHs(0)
    sulfur_atom.SetNoImplicit(True)

    # Check if the molecule contains the CoA substructure
    if not mol.HasSubstructMatch(coa_mol):
        return False, "Coenzyme A moiety not found in molecule"

    # Define thioester pattern with sulfur connected to carbonyl group
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")

    # Find thioester linkage in the molecule
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Verify that sulfur in thioester is part of CoA moiety
    coa_match = mol.GetSubstructMatch(coa_mol)
    if not coa_match:
        return False, "Coenzyme A moiety not properly matched"

    # Map CoA sulfur atom to molecule sulfur atom
    coa_sulfur_idx = sulfur_idx
    mol_sulfur_idx = coa_match[coa_sulfur_idx]

    # Check that sulfur atom in thioester is the same as in CoA moiety
    thioester_sulfur_idxs = [match[1] for match in thioester_matches]
    if mol_sulfur_idx not in thioester_sulfur_idxs:
        return False, "Thioester linkage not formed with CoA sulfur atom"

    # Determine fatty acyl chain length
    # Start from carbonyl carbon of thioester linkage
    for match in thioester_matches:
        if match[1] == mol_sulfur_idx:
            carbonyl_c_idx = match[0]
            break
    else:
        return False, "Carbonyl carbon not found for thioester linkage"

    # Traverse fatty acyl chain starting from carbonyl carbon
    fatty_acyl_chain = set()
    atoms_to_visit = [mol.GetAtomWithIdx(carbonyl_c_idx)]
    visited_atoms = set()
    while atoms_to_visit:
        atom = atoms_to_visit.pop()
        atom_idx = atom.GetIdx()
        if atom_idx in visited_atoms or atom_idx in coa_match:
            continue
        visited_atoms.add(atom_idx)
        fatty_acyl_chain.add(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited_atoms and neighbor_idx not in coa_match:
                atoms_to_visit.append(neighbor)

    # Count number of carbons in fatty acyl chain (excluding carbonyl carbon)
    num_carbons = sum(1 for idx in fatty_acyl_chain if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if num_carbons < 4:
        return False, f"Fatty acyl chain too short ({num_carbons} carbons)"

    return True, "Molecule is a fatty acyl-CoA with appropriate Coenzyme A structure"