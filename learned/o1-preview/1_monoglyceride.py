"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is a monoglyceride in which the acyl substituent is located at position 1,
    and positions 2 and 3 have free hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove salts and separate fragments (keep the largest)
    mol_frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(mol_frags) > 1:
        mol = max(mol_frags, default=mol, key=lambda m: m.GetNumAtoms())

    # Define the ester functional group
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected exactly 1 ester group, found {len(ester_matches)}"

    # Define the glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    glycerol_matches = mol.GetSubstructMatch(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone with free hydroxyls at positions 2 and 3 found"

    # Map glycerol atoms to positions
    glycerol_atoms = [mol.GetAtomWithIdx(idx) for idx in glycerol_matches]

    # Identify the esterified hydroxyl at position 1
    position_1 = None
    position_2 = None
    position_3 = None
    for atom in glycerol_atoms:
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 6 and neighbor.IsInRing() == False:
                    # Check if this oxygen is involved in an ester linkage
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        attached_carbon = neighbor
                        attached_carbon_neighbors = attached_carbon.GetNeighbors()
                        for acn in attached_carbon_neighbors:
                            if acn.GetAtomicNum() == 6 and acn != atom:
                                # Found a carbon chain attached to this oxygen via carbon
                                ester_bond = mol.GetBondBetweenAtoms(attached_carbon.GetIdx(), acn.GetIdx())
                                if ester_bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                    # This oxygen is esterified
                                    position_1 = atom
                                else:
                                    position_2 = atom
                            else:
                                position_3 = atom

    if not position_1:
        return False, "No esterified hydroxyl at position 1 found"

    if not position_2 or not position_3:
        return False, "Positions 2 and/or 3 do not have free hydroxyl groups"

    # Check for absence of other functional groups (e.g., phosphates)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Molecule contains a phosphate group"

    # Check for other ester or amide groups
    other_ester_pattern = Chem.MolFromSmarts("[$([CX3](=O)[OX2H1]),$([CX3](=O)[OX1-])]")
    other_esters = mol.GetSubstructMatches(other_ester_pattern)
    if len(other_esters) > 1:
        return False, "Molecule contains additional ester or amide groups"

    return True, "Molecule matches 1-monoglyceride structure with acyl group at position 1 and free hydroxyls at positions 2 and 3"