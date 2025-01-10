"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol is a glycerophosphoinositol having one phosphatidyl group
    esterified to one of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Normalize molecule (add hydrogens)
    mol = Chem.AddHs(mol)

    # Generate inositol ring molecule (without stereochemistry)
    inositol_smiles = 'OC1C(O)C(O)C(O)C(O)C1O'
    inositol_mol = Chem.MolFromSmiles(inositol_smiles)

    # Check for inositol substructure
    if not mol.HasSubstructMatch(inositol_mol):
        return False, "Inositol ring not found"

    # Check for phosphate group connected to inositol
    # Phosphate group connected via oxygen to inositol
    phosphate_pattern = Chem.MolFromSmarts('O[P](=O)(O)O')
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)

    # Find phosphate groups connected to inositol oxygen
    phosphate_bonded_to_inositol = False
    for match in phosphate_matches:
        phosphate_idx = match[1]  # Index of phosphorus atom
        phosphate_atom = mol.GetAtomWithIdx(phosphate_idx)
        neighbors = phosphate_atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.HasSubstructMatch(inositol_mol):
                phosphate_bonded_to_inositol = True
                break
        if phosphate_bonded_to_inositol:
            break
    if not phosphate_bonded_to_inositol:
        return False, "Phosphate group connected to inositol not found"

    # Check for glycerol backbone connected to phosphate
    glycerol_smiles = 'OCC(O)CO'
    glycerol_mol = Chem.MolFromSmiles(glycerol_smiles)
    if not mol.HasSubstructMatch(glycerol_mol):
        return False, "Glycerol backbone not found"

    # Check if glycerol is connected to phosphate
    glycerol_atoms = mol.GetSubstructMatch(glycerol_mol)
    phosphate_atoms = mol.GetSubstructMatch(phosphate_pattern)
    connected = False
    for gid in glycerol_atoms:
        glycerol_atom = mol.GetAtomWithIdx(gid)
        for pid in phosphate_atoms:
            phosphate_atom = mol.GetAtomWithIdx(pid)
            bond = mol.GetBondBetweenAtoms(gid, pid)
            if bond:
                connected = True
                break
        if connected:
            break
    if not connected:
        return False, "Glycerol not connected to phosphate group"

    # Check for ester linkages at sn-1 and sn-2 positions
    ester_pattern = Chem.MolFromSmarts('C(=O)O[C;H]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, "Less than 2 ester bonds found at sn-1 and sn-2 positions"

    return True, "Molecule matches phosphatidylinositol structure"