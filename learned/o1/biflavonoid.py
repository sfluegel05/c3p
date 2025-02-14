"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: biflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is a flavonoid oligomer obtained by the oxidative coupling of at least two flavonoid units (aryl-substituted benzopyran rings or derivatives), resulting in two ring systems being joined together by a single atom or bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define flavonoid core SMARTS pattern (chromone moiety)
    flavonoid_core_smarts = 'c1cc2oc(=O)cc2cc1'  # chromone-like core
    flavonoid_core = Chem.MolFromSmarts(flavonoid_core_smarts)
    if flavonoid_core is None:
        return None, "Invalid flavonoid core SMARTS pattern"

    # Find matches for the flavonoid core
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_core)
    num_flavonoid_cores = len(flavonoid_matches)

    if num_flavonoid_cores < 2:
        return False, f"Found {num_flavonoid_cores} flavonoid core(s), need at least 2"

    # Create sets of atoms for each flavonoid core
    flavonoid_core_atoms = [set(match) for match in flavonoid_matches]

    # Find pairs of flavonoid cores that are connected
    connected_pairs = []
    for i in range(num_flavonoid_cores):
        for j in range(i+1, num_flavonoid_cores):
            # Check if there is a bond between any atom of core i and any atom of core j
            core_i_atoms = flavonoid_core_atoms[i]
            core_j_atoms = flavonoid_core_atoms[j]
            for atom_i in core_i_atoms:
                for atom_j in core_j_atoms:
                    bond = mol.GetBondBetweenAtoms(atom_i, atom_j)
                    if bond is not None:
                        # Found a bond connecting the two cores
                        connected_pairs.append((i, j))
                        break
                if (i, j) in connected_pairs:
                    break

    if len(connected_pairs) == 0:
        return False, "Flavonoid cores are not connected"

    # Check if the connection is via a single atom or bond
    for (i, j) in connected_pairs:
        num_bonds_between_cores = 0
        core_i_atoms = flavonoid_core_atoms[i]
        core_j_atoms = flavonoid_core_atoms[j]
        for atom_i in core_i_atoms:
            for atom_j in core_j_atoms:
                bond = mol.GetBondBetweenAtoms(atom_i, atom_j)
                if bond is not None:
                    num_bonds_between_cores += 1
        if num_bonds_between_cores == 1:
            return True, "Contains at least two flavonoid cores connected via a single bond"
        else:
            return False, f"Flavonoid cores are connected via {num_bonds_between_cores} bonds, expected 1"

    return False, "Could not determine connectivity between flavonoid cores"


__metadata__ = {   'chemical_class': {   'name': 'biflavonoid',
                              'definition': 'A flavonoid oligomer that is obtained by the oxidative coupling of at least two units of aryl-substituted benzopyran rings or its substituted derivatives, resulting in the two ring systems being joined together by a single atom or bond.'},
        'config': {   'llm_model_name': 'lbl/claude-sonnet'},
        'success': True}