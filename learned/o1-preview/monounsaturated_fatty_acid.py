"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: monounsaturated fatty acid
"""
from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    A monounsaturated fatty acid (MUFA) is a fatty acid with exactly one double or triple bond
    in the fatty acid chain and singly bonded carbon atoms in the rest of the chain.
    The fatty acid chain is defined as the longest unbranched carbon chain starting from a
    carboxyl carbon (including esterified forms).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid or ester group (allowing for esterified fatty acids)
    # SMARTS pattern for carboxylic acid: C(=O)[OH]
    # SMARTS pattern for ester: C(=O)O
    functional_group = Chem.MolFromSmarts('C(=O)[O;H1,H0,-1]')  # Carboxylic acid or deprotonated form
    ester_group = Chem.MolFromSmarts('C(=O)O[CX4]')  # Ester linkage
    groups = mol.GetSubstructMatches(functional_group) + mol.GetSubstructMatches(ester_group)
    
    if not groups:
        return False, "No carboxylic acid or ester group found"

    # Attempt to find a fatty acid chain starting from each functional group
    for group in groups:
        carbonyl_c_idx = group[0]  # Carbonyl carbon index
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)

        # Find the alpha carbon (carbon attached to carbonyl carbon)
        alpha_c = None
        for nbr in carbonyl_c.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in group:
                alpha_c = nbr
                break
        if alpha_c is None:
            continue  # No alpha carbon found; try next match

        visited = set()
        double_bonds = 0
        branching = False

        # Recursive function to traverse the fatty acid chain
        def traverse_chain(atom, prev_atom):
            nonlocal double_bonds, branching
            atom_idx = atom.GetIdx()
            if atom_idx in visited:
                return
            visited.add(atom_idx)

            # Count the number of carbon neighbors (excluding the previous atom)
            carbon_nbrs = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom.GetIdx()]
            if len(carbon_nbrs) > 1:
                branching = True  # More than two carbon neighbors indicates branching

            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetIdx() == prev_atom.GetIdx():
                    continue  # Skip the bond we just came from
                if nbr.GetAtomicNum() != 6:
                    continue  # Skip non-carbon atoms

                bond_type = bond.GetBondType()
                # Count double or triple bonds between carbons
                if bond_type in (Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE):
                    double_bonds += 1
                traverse_chain(nbr, atom)

        # Start traversing from the alpha carbon
        traverse_chain(alpha_c, carbonyl_c)

        # Check for branching and double bond count
        if branching:
            continue  # Not an unbranched chain; try next match
        if double_bonds != 1:
            continue  # Does not have exactly one double/triple bond; try next match

        return True, "Molecule is a monounsaturated fatty acid"

    return False, "Does not match criteria for a monounsaturated fatty acid"