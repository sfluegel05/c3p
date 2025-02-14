"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolops

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    An unsaturated fatty acyl-CoA results from the formal condensation of the thiol group
    of coenzyme A with the carboxy group of any unsaturated fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coenzyme A core SMARTS pattern (flexible to accommodate variations)
    coa_smarts = "NC(=O)CCNC(=O)[C@H](O)[C@H](C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](CO[P](=O)(O)O)[C@@H](O)[C@H]1O"
    coa_mol = Chem.MolFromSmarts(coa_smarts)
    if coa_mol is None:
        return False, "Error generating CoA pattern from SMARTS"

    # Find CoA substructure matches
    coa_matches = mol.GetSubstructMatches(coa_mol, useChirality=False)
    if not coa_matches:
        return False, "No Coenzyme A moiety found"

    # Define thioester linkage pattern (C(=O)-S)
    thioester_smarts = "C(=O)S"
    thioester_mol = Chem.MolFromSmarts(thioester_smarts)
    thioester_matches = mol.GetSubstructMatches(thioester_mol)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Identify the fatty acyl chain connected via thioester bond to CoA
    for thioester_match in thioester_matches:
        carbonyl_c_idx = thioester_match[0]  # Carbonyl carbon index
        sulfur_idx = thioester_match[2]      # Sulfur atom index

        # Check if sulfur is connected to CoA
        is_sulfur_in_coa = False
        for coa_match in coa_matches:
            if sulfur_idx in coa_match:
                is_sulfur_in_coa = True
                break
        if not is_sulfur_in_coa:
            continue  # Sulfur not in CoA, skip to next thioester linkage

        # Extract fatty acyl chain (atoms connected to carbonyl carbon, excluding CoA)
        fatty_acyl_atoms = set()
        atoms_to_explore = [carbonyl_c_idx]
        while atoms_to_explore:
            atom_idx = atoms_to_explore.pop()
            if atom_idx in fatty_acyl_atoms:
                continue
            fatty_acyl_atoms.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Skip sulfur atom and CoA atoms
                if neighbor_idx == sulfur_idx:
                    continue
                is_neighbor_in_coa = any(neighbor_idx in coa_match for coa_match in coa_matches)
                if is_neighbor_in_coa:
                    continue
                if neighbor_idx not in fatty_acyl_atoms:
                    atoms_to_explore.append(neighbor_idx)

        # Create sub-molecule of the fatty acyl chain
        fatty_acyl_submol = Chem.PathToSubmol(mol, list(fatty_acyl_atoms))
        if fatty_acyl_submol is None:
            continue

        # Check for unsaturation (presence of carbon-carbon double bonds)
        has_double_bond = False
        for bond in fatty_acyl_submol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                # Check if both atoms are carbons
                if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                    has_double_bond = True
                    break

        if has_double_bond:
            return True, (
                "Contains CoA moiety with unsaturated fatty acyl chain "
                "linked via thioester bond"
            )
        else:
            return False, "Fatty acyl chain is saturated (no carbon-carbon double bonds)"

    return False, "No thioester linkage connected to CoA was found"