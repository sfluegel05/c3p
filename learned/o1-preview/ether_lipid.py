"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: CHEBI:12307 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid is similar in structure to a glycerolipid but in which one or more of the carbon atoms on glycerol
    is bonded to an alkyl chain via an ether linkage, as opposed to the usual ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (three carbons connected in sequence)
    glycerol_pattern = Chem.MolFromSmarts("C-C-C")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # For each glycerol backbone found
    for match in matches:
        c1_idx, c2_idx, c3_idx = match
        glycerol_carbon_indices = set([c1_idx, c2_idx, c3_idx])

        # Check for ether linkage on glycerol carbons
        ether_found = False
        for c_idx in glycerol_carbon_indices:
            atom = mol.GetAtomWithIdx(c_idx)
            # Iterate over neighbors to find oxygen connected via single bond
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(c_idx, neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        oxygen_atom = neighbor
                        # Check if oxygen is connected to another carbon not in glycerol backbone
                        for o_neighbor in oxygen_atom.GetNeighbors():
                            if o_neighbor.GetIdx() != c_idx and o_neighbor.GetAtomicNum() == 6:
                                if o_neighbor.GetIdx() not in glycerol_carbon_indices:
                                    o_bond = mol.GetBondBetweenAtoms(oxygen_atom.GetIdx(), o_neighbor.GetIdx())
                                    if o_bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                        # Found ether linkage
                                        ether_found = True
                                        break
                        if ether_found:
                            break
            if ether_found:
                break

        if ether_found:
            return True, "Contains glycerol backbone with ether linkage"

    return False, "No ether linkage found on glycerol backbone"