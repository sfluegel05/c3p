"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for the carboxylic acid group pattern
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylic acid group found"

    # Start from the carboxylic group and trace the carbon chain
    for match in mol.GetSubstructMatches(carboxylate_pattern):
        carboxyl_carbon = mol.GetAtomWithIdx(match[0])

        # Traverse the longest straight carbon chain starting near carboxyl group
        visited_atoms = set()
        stack = [(carboxyl_carbon, None)]  # atom and previous carbon

        while stack:
            atom, prev_carbon = stack.pop()
            visited_atoms.add(atom.GetIdx())

            # Check saturation and branching:
            if atom.GetAtomicNum() == 6:  # Carbon
                # non-terminal carbons past the first should have exactly two hydrogen removals
                single_bond_partners = [
                    b.GetOtherAtom(atom) for b in atom.GetBonds() 
                    if b.GetBondType() == Chem.BondType.SINGLE and b.GetOtherAtom(atom).GetIdx() != prev_carbon
                ]

                if len(single_bond_partners) > 2:
                    return False, "Carbon chain is branched"

                for bond in atom.GetBonds():
                    # Verify it's part of the main chain
                    if bond.GetBondType() not in (Chem.BondType.SINGLE,):
                        return False, "Carbon chain contains double or triple bonds"

                # Continue visiting neighbors
                for partner in single_bond_partners:
                    if partner.GetIdx() not in visited_atoms:
                        stack.append((partner, atom.GetIdx()))

    return True, "Molecule is a straight-chain saturated fatty acid"

# Examples
example_smiles = "CCCCCCC(O)=O"  # Heptanoic acid
result, reason = is_straight_chain_saturated_fatty_acid(example_smiles)
print(f"{example_smiles}: {result} ({reason})")