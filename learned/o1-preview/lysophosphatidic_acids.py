"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem

def is_lysophosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    Lysophosphatidic acids are monoacylglycerol phosphates obtained by hydrolytic removal 
    of one of the two acyl groups of any phosphatidic acid or derivatives therein.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find phosphorus atoms (phosphate groups)
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_atoms:
        return False, "No phosphate group found"

    # Assume only one phosphate group for LPA
    phosphate_atom = phosphorus_atoms[0]

    # Find the oxygen atom connecting phosphate to glycerol backbone
    glycerol_phosphate_oxygen = None
    for neighbor in phosphate_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8:
            for nbr in neighbor.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != phosphate_atom.GetIdx():
                    glycerol_phosphate_oxygen = neighbor
                    break
            if glycerol_phosphate_oxygen:
                break

    if not glycerol_phosphate_oxygen:
        return False, "No phosphate ester linkage found"

    # Get the first carbon of glycerol backbone connected to phosphate
    glycerol_carbon1 = None
    for atom in glycerol_phosphate_oxygen.GetNeighbors():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() != phosphate_atom.GetIdx():
            glycerol_carbon1 = atom
            break

    if not glycerol_carbon1:
        return False, "Cannot identify glycerol backbone starting point"

    # Trace the glycerol backbone (three carbons connected in sequence)
    glycerol_carbons = [glycerol_carbon1]
    current_atom = glycerol_carbon1
    prev_atom = glycerol_phosphate_oxygen
    while len(glycerol_carbons) < 3:
        neighbors = [atom for atom in current_atom.GetNeighbors() if atom.GetAtomicNum() == 6 and atom.GetIdx() != prev_atom.GetIdx()]
        if not neighbors:
            return False, "Incomplete glycerol backbone"
        next_atom = neighbors[0]
        glycerol_carbons.append(next_atom)
        prev_atom = current_atom
        current_atom = next_atom

    glycerol_carbon_indices = [atom.GetIdx() for atom in glycerol_carbons]

    # Count ester bonds connected to glycerol carbons (acyl chains)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    acyl_ester_count = 0
    for match in ester_matches:
        ester_oxygen_idx = match[1]
        connected_carbon_idx = match[2]
        if connected_carbon_idx in glycerol_carbon_indices:
            acyl_ester_count += 1

    if acyl_ester_count != 1:
        return False, f"Expected 1 acyl chain attached via ester bond to glycerol backbone, found {acyl_ester_count}"

    # Check that remaining glycerol carbons have correct substituents
    # One should be connected to phosphate via oxygen (already confirmed)
    # One should have a hydroxyl group
    hydroxyl_count = 0
    for carbon_atom in glycerol_carbons:
        carbon_idx = carbon_atom.GetIdx()
        if carbon_idx in glycerol_carbon_indices:
            # Check for hydroxyl groups (O connected to H)
            for neighbor in carbon_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    # Ensure oxygen is connected to hydrogen
                    has_hydrogen = False
                    for nbr in neighbor.GetNeighbors():
                        if nbr.GetAtomicNum() == 1:
                            has_hydrogen = True
                            break
                    if has_hydrogen:
                        hydroxyl_count += 1
                        break

    if hydroxyl_count != 1:
        return False, f"Expected 1 free hydroxyl group on glycerol backbone, found {hydroxyl_count}"

    return True, "Molecule is a lysophosphatidic acid with glycerol phosphate backbone and one acyl chain"