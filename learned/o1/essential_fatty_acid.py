"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: CHEBI:18035 essential fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    An essential fatty acid is a polyunsaturated fatty acid that cannot be synthesized
    by the human body and must be obtained from the diet.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all carboxylic acid groups (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    if not carboxylic_acid_matches:
        return False, "No carboxylic acid groups found"

    essential = False
    reasons = []

    for match in carboxylic_acid_matches:
        # The carboxylic acid carbon atom
        carboxyl_carbon_idx = match[0]
        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

        # Start traversal from the carboxyl carbon to find the longest chain
        visited = set()
        chain_atoms = []
        double_bond_positions = []
        branches = False

        def dfs(atom, position):
            nonlocal branches
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() == 6 and not atom.IsInRing():
                chain_atoms.append(atom.GetIdx())
                for bond in atom.GetBonds():
                    neighbor = bond.GetOtherAtom(atom)
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx == carboxyl_carbon_idx:
                        continue  # Do not go back to carboxyl carbon
                    if neighbor_idx not in visited:
                        if neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing():
                            if len(atom.GetNeighbors()) > 2:
                                branches = True  # Detected a branch
                            # Record double bond positions (delta numbering from carboxyl carbon)
                            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                                double_bond_positions.append(position)
                            dfs(neighbor, position + 1)

        # Start DFS from the carbon connected to carboxyl carbon (if any)
        for bond in carboxyl_carbon.GetBonds():
            neighbor = bond.GetOtherAtom(carboxyl_carbon)
            if neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing():
                if neighbor.GetIdx() not in visited:
                    if len(carboxyl_carbon.GetNeighbors()) > 2:
                        branches = True
                    if bond.GetBondType() == rdchem.BondType.DOUBLE:
                        double_bond_positions.append(1)
                    dfs(neighbor, 1)

        num_carbons = len(chain_atoms) + 1  # Include the carboxyl carbon
        num_double_bonds = len(double_bond_positions)

        if branches:
            reasons.append(f"Chain has branches, not a linear fatty acid")
            continue  # Skip branched chains

        # Essential fatty acid criteria
        if num_carbons >= 16 and num_double_bonds >= 2:
            # Check for double bonds beyond delta-9 position
            beyond_delta_9 = [pos for pos in double_bond_positions if pos > 9]
            if beyond_delta_9:
                essential = True
                reasons.append(
                    f"Found fatty acid chain with {num_carbons} carbons and "
                    f"double bonds at positions {double_bond_positions} (delta numbering)"
                )
            else:
                reasons.append(
                    f"Fatty acid chain with {num_carbons} carbons but no double bonds beyond delta-9"
                )
        else:
            reasons.append(
                f"Fatty acid chain with {num_carbons} carbons and {num_double_bonds} double bonds does not meet criteria"
            )

    if essential:
        return True, "Molecule contains essential fatty acid chain(s): " + "; ".join(reasons)
    else:
        return False, "No essential fatty acid chains found: " + "; ".join(reasons)