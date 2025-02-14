"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: CHEBI:18035 essential fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdMolDescriptors

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

    # Define a SMARTS pattern for fatty acyl chains (aliphatic chain with optional ester linkage)
    fatty_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[O,Z][C;!$(C=O);!$(C=O[O,N])].[CX4;R0][CX4;R0]')
    # The pattern matches acyl groups connected via ester or amide linkages to aliphatic chains

    # Find all matches of the fatty acid pattern
    matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not matches:
        return False, "No fatty acid chains found"

    essential = False
    reasons = []
    # Process each fatty acid chain
    for match in matches:
        carbon_chain = set()
        double_bond_positions = []
        visited = set()
        queue = []

        # Get the atom index of the carbonyl carbon
        carbonyl_carbon_idx = match[0]
        # Get the atom index of the alpha carbon (next carbon in the chain)
        alpha_carbon_idx = match[2]

        # Start traversal from the alpha carbon
        queue.append((alpha_carbon_idx, 1))  # (atom_idx, position in chain)

        while queue:
            atom_idx, position = queue.pop(0)
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)

            if atom.GetAtomicNum() == 6 and not atom.IsInRing():
                carbon_chain.add(atom_idx)
                # Check for double bonds
                for bond in atom.GetBonds():
                    neighbor = bond.GetOtherAtom(atom)
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited:
                        if bond.GetBondType() == rdchem.BondType.DOUBLE and neighbor.GetAtomicNum() == 6:
                            # Record position of the double bond relative to alpha carbon
                            double_bond_positions.append(position)
                        if neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing():
                            queue.append((neighbor_idx, position + 1))

        num_carbons = len(carbon_chain) + 1  # Include the carbonyl carbon
        num_double_bonds = len(set(double_bond_positions))

        # Check for essential fatty acid criteria
        if num_carbons >= 16 and num_double_bonds >= 2:
            # Check if any double bonds are beyond delta-9 position
            beyond_delta_9 = any(pos > 9 for pos in double_bond_positions)
            if beyond_delta_9:
                essential = True
                reasons.append(
                    f"Found fatty acid chain with {num_carbons} carbons and "
                    f"{num_double_bonds} double bonds beyond delta-9 position"
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