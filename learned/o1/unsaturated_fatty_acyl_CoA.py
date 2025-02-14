"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdmolops import GetShortestPath

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

    # Coenzyme A SMILES string from PubChem (CID 977)
    coa_smiles = (
        "C[C@@H](O)C(=O)NCCSCCNC(=O)CCNC(=O)"
        "[C@H](O)[C@H](C)(C)COP(=O)(O)OP(=O)(O)"
        "OC[C@H]1O[C@H]([C@@H](O)[C@H]1O)"
        "N1C=NC2=C1N=CN=C2N"
    )
    coa_mol = Chem.MolFromSmiles(coa_smiles)
    if coa_mol is None:
        return False, "Error generating CoA molecule from SMILES"

    # Perform substructure search for CoA in the molecule without chirality
    if not mol.HasSubstructMatch(coa_mol, useChirality=False):
        return False, "No Coenzyme A moiety found"

    # Define thioester linkage pattern (C(=O)-S-C)
    thioester_smarts = "C(=O)S"
    thioester_mol = Chem.MolFromSmarts(thioester_smarts)
    if thioester_mol is None:
        return False, "Error parsing thioester SMARTS pattern"

    # Find thioester linkages in the molecule
    thioester_matches = mol.GetSubstructMatches(thioester_mol)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # For each thioester linkage, check if sulfur is connected to CoA
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # Index of the carbonyl carbon
        sulfur_idx = match[2]      # Index of the sulfur atom

        # Get indices of CoA atoms in the molecule
        coa_matches = mol.GetSubstructMatches(coa_mol, useChirality=False)
        coa_atom_indices = set()
        for cm in coa_matches:
            coa_atom_indices.update(cm)
        if not coa_atom_indices:
            continue  # Should not happen as CoA is present

        # Check if sulfur atom is connected to CoA
        paths_to_coa = [
            GetShortestPath(mol, sulfur_idx, coa_idx) for coa_idx in coa_atom_indices
        ]
        paths_to_coa = [path for path in paths_to_coa if path]
        if not paths_to_coa:
            continue  # Sulfur not connected to CoA, check next thioester linkage

        # Extract fatty acyl chain connected to the carbonyl carbon
        fatty_acyl_atoms = set()
        visited = set()
        to_visit = [carbonyl_c_idx]
        while to_visit:
            current_idx = to_visit.pop()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            fatty_acyl_atoms.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Exclude sulfur and CoA atoms
                if neighbor_idx == sulfur_idx or neighbor_idx in coa_atom_indices:
                    continue
                if neighbor_idx not in visited:
                    to_visit.append(neighbor_idx)

        # Create a sub-molecule of the fatty acyl chain
        fatty_acyl_submol = Chem.PathToSubmol(mol, fatty_acyl_atoms)

        # Check for unsaturation (presence of C=C bonds)
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