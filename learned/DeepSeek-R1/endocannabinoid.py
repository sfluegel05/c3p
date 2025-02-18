"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: endocannabinoid (CHEBI:58997)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids are characterized by ethanolamide groups or glycerol esters/ethers with long fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Function to calculate chain length iteratively
    def get_chain_length(start_atom):
        length = 0
        current_atom = start_atom
        visited = set()
        while current_atom and current_atom.GetAtomicNum() == 6 and current_atom.GetIdx() not in visited:
            visited.add(current_atom.GetIdx())
            length += 1
            # Get next carbon in chain (non-aromatic, non-ring, single/double bond)
            next_atoms = []
            for neighbor in current_atom.GetNeighbors():
                if (neighbor.GetAtomicNum() == 6 and
                    not neighbor.GetIsAromatic() and
                    not neighbor.IsInRing() and
                    neighbor.GetIdx() not in visited):
                    bond = mol.GetBondBetweenAtoms(current_atom.GetIdx(), neighbor.GetIdx())
                    if bond and bond.GetBondType() in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE):
                        next_atoms.append(neighbor)
            if len(next_atoms) != 1:
                break  # Branch or end
            current_atom = next_atoms[0]
        return length

    # Check for ethanolamide structure (NCCO connected via amide)
    ethanolamide_pattern = Chem.MolFromSmarts("[NH]CCO")
    ethanolamide_matches = mol.GetSubstructMatches(ethanolamide_pattern)
    for match in ethanolamide_matches:
        n_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in neighbor.GetBonds()):
                carbonyl_c = neighbor
                r_group = [a for a in carbonyl_c.GetNeighbors() if a != n_atom]
                if r_group:
                    chain_length = get_chain_length(r_group[0])
                    if chain_length >= 12:
                        return True, "Ethanolamide with long fatty acid chain (>=12 carbons)"

    # Check for glycerol structure with ester/ether and long chain
    glycerol_backbone = Chem.MolFromSmarts("C-C-C")
    backbone_matches = mol.GetSubstructMatches(glycerol_backbone)
    for match in backbone_matches:
        c1, c2, c3 = [mol.GetAtomWithIdx(idx) for idx in match]
        backbone_atoms = {c1.GetIdx(), c2.GetIdx(), c3.GetIdx()}
        hydroxyl_count = 0
        ether_ester_oxygen = None

        # Check substituents on each carbon in the backbone
        for carbon in [c1, c2, c3]:
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    # Check for hydroxyl
                    if neighbor.GetTotalNumHs() >= 1 and neighbor.GetDegree() == 1:
                        hydroxyl_count += 1
                    else:
                        # Check for ester or ether
                        ester = False
                        for nbr in neighbor.GetNeighbors():
                            if nbr.GetAtomicNum() == 6 and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in nbr.GetBonds()):
                                ester = True
                                break
                        if ester:
                            # Ester found, get carbonyl carbon
                            for nbr in neighbor.GetNeighbors():
                                if any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in nbr.GetBonds()):
                                    r_group = [a for a in nbr.GetNeighbors() if a != neighbor]
                                    if r_group:
                                        chain_length = get_chain_length(r_group[0])
                                        if chain_length >= 12:
                                            return True, "Glycerol ester with long chain (>=12 carbons)"
                        else:
                            # Check for ether
                            if neighbor.GetDegree() == 2:
                                other_carbons = [n for n in neighbor.GetNeighbors() if n.GetIdx() not in backbone_atoms]
                                if other_carbons:
                                    chain_length = get_chain_length(other_carbons[0])
                                    if chain_length >= 12:
                                        return True, "Glycerol ether with long chain (>=12 carbons)"

        # Check if two hydroxyls and one ether/ester found
        if hydroxyl_count >= 2 and ether_ester_oxygen:
            return True, "Glycerol backbone with two hydroxyls and one long-chain ether/ester"

    return False, "Does not match endocannabinoid structural criteria"