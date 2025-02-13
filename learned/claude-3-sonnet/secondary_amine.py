"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:35612 secondary amine
A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of nitrogen atom
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atom found"

    # Check for nitrogen atoms bonded to two alkyl/aryl groups
    is_secondary_amine = False
    for n_atom in nitrogen_atoms:
        n_neighbors = [neighbor for neighbor in n_atom.GetNeighbors() if neighbor.GetAtomicNum() == 6]
        if len(n_neighbors) >= 2:
            # Identify alkyl/aryl groups
            alkyl_aryl_groups = []
            for neighbor in n_neighbors:
                alkyl_aryl_group = [neighbor]
                visited = set()
                to_visit = [neighbor]
                while to_visit:
                    cur_atom = to_visit.pop(0)
                    if cur_atom.GetAtomicNum() == 6:
                        alkyl_aryl_group.append(cur_atom)
                        visited.add(cur_atom)
                        for neighbor in cur_atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 6 and neighbor not in visited:
                                to_visit.append(neighbor)
                alkyl_aryl_groups.append(alkyl_aryl_group)

            # Check for at least two separate alkyl/aryl groups
            if len(alkyl_aryl_groups) >= 2:
                is_secondary_amine = True
                break

    if not is_secondary_amine:
        return False, "Nitrogen atom not bonded to at least two alkyl/aryl groups"

    # Check for common functional groups/heteroatoms
    heteroatoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6, 7]]
    if heteroatoms:
        functional_groups = []
        for atom in heteroatoms:
            if atom.GetAtomicNum() == 8:  # Oxygen
                functional_groups.append("hydroxy" if sum(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors()) else "carbonyl")
            elif atom.GetAtomicNum() == 16:  # Sulfur
                functional_groups.append("thiol" if sum(neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors()) else "sulfur")
            else:
                functional_groups.append(f"heteroatom {atom.GetSymbol()}")
        reason = f"Contains a nitrogen atom bonded to at least two alkyl/aryl groups, with additional functional group(s): {', '.join(functional_groups)}"
    else:
        reason = "Contains a nitrogen atom bonded to at least two alkyl/aryl groups"

    return True, reason