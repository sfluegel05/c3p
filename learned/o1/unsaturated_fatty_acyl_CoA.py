"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Define Coenzyme A (CoA) substructure using a simplified SMARTS pattern
    # This pattern captures the phosphoadenosine diphosphate part
    coA_smarts = (
        "O=P(O)(O)OC[C@H]1O[C@H](CO[P](O)(O)=O)[C@@H](O)[C@H]1O"
        "n1cnc2c(N)ncnc12"
    )
    coA_mol = Chem.MolFromSmarts(coA_smarts)
    if coA_mol is None:
        return False, "Error parsing CoA SMARTS pattern"
    
    if not mol.HasSubstructMatch(coA_mol):
        return False, "No Coenzyme A moiety found"

    # Define thioester linkage pattern (C(=O)-S)
    thioester_smarts = "C(=O)S"
    thioester_mol = Chem.MolFromSmarts(thioester_smarts)
    if thioester_mol is None:
        return False, "Error parsing thioester SMARTS pattern"

    # Find thioester linkages
    thioester_matches = mol.GetSubstructMatches(thioester_mol)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # For each thioester linkage, check if it's connected to CoA
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon index
        sulfur_idx = match[1]      # Sulfur atom index

        # Use shortest path to determine if sulfur atom is connected to CoA
        coa_matches = mol.GetSubstructMatch(coA_mol)
        coa_atom_indices = set(coa_matches)
        if not coa_atom_indices:
            continue  # No CoA matched, but should not happen due to earlier check

        paths = Chem.rdmolops.GetShortestPath(mol, sulfur_idx, list(coa_atom_indices)[0])
        if not paths:
            continue  # Sulfur not connected to CoA

        # Extract the fatty acyl chain (atoms connected to carbonyl carbon excluding CoA)
        # We will perform a BFS starting from the carbonyl carbon and stopping at sulfur and CoA
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
                if neighbor_idx == sulfur_idx:
                    continue  # Do not cross thioester bond to sulfur
                if neighbor_idx in coa_atom_indices:
                    continue  # Do not include CoA atoms
                if neighbor_idx not in visited:
                    to_visit.append(neighbor_idx)

        # Create a sub-molecule of the fatty acyl chain
        fatty_acyl_submol = Chem.PathToSubmol(mol, fatty_acyl_atoms)

        # Check for unsaturation (C=C bonds)
        has_double_bond = False
        for bond in fatty_acyl_submol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                    has_double_bond = True
                    break

        if has_double_bond:
            return True, "Contains CoA moiety with unsaturated fatty acyl chain linked via thioester bond"
        else:
            return False, "Fatty acyl chain is saturated (no carbon-carbon double bonds)"

    return False, "No thioester linkage connected to CoA was found"