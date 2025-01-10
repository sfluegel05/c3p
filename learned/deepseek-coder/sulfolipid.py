"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: CHEBI:63505 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is a compound containing a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfonic acid group (-SO3H)
    sulfonic_acid_pattern = Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])([OX2H])")
    if not mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "No sulfonic acid group found"

    # Check for a carbon-sulfur bond (C-S) connected to the sulfonic acid group
    carbon_sulfur_bond_pattern = Chem.MolFromSmarts("[CX4]-[SX4](=[OX1])(=[OX1])([OX2H])")
    if not mol.HasSubstructMatch(carbon_sulfur_bond_pattern):
        return False, "No carbon-sulfur bond connected to sulfonic acid group"

    # Check for lipid-like structure (long hydrocarbon chains or ester/amide bonds)
    # We can approximate this by checking for a minimum number of carbons and a certain ratio of carbons to other atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)

    # A lipid typically has a high carbon-to-oxygen ratio
    if c_count < 20 or c_count / (o_count + s_count) < 2:
        return False, "Molecule does not have a lipid-like structure (insufficient carbon content)"

    # Check for long hydrocarbon chains (at least 8 carbons in a row)
    long_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        # If no long chain, check for ester or amide bonds
        ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]")
        amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3][CX4]")
        if not (mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(amide_pattern)):
            return False, "No long hydrocarbon chain or ester/amide bond found"

    # Check molecular weight - sulfolipids typically have a higher molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for sulfolipid"

    # Ensure the sulfonic acid group is directly connected to a lipid-like structure
    # by checking that the carbon connected to sulfur is part of a long chain or ester/amide bond
    sulfur_atom = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16 and atom.GetDegree() == 4:  # Sulfur in sulfonic acid
            sulfur_atom = atom
            break

    if sulfur_atom is None:
        return False, "No sulfonic acid group found"

    carbon_neighbor = None
    for neighbor in sulfur_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Carbon
            carbon_neighbor = neighbor
            break

    if carbon_neighbor is None:
        return False, "No carbon-sulfur bond found"

    # Check if the carbon connected to sulfur is part of a long chain or ester/amide bond
    long_chain_carbon_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    ester_carbon_pattern = Chem.MolFromSmarts("[CX4]-[CX3](=[OX1])[OX2][CX4]")
    amide_carbon_pattern = Chem.MolFromSmarts("[CX4]-[CX3](=[OX1])[NX3][CX4]")
    if not (mol.HasSubstructMatch(long_chain_carbon_pattern) or 
            mol.HasSubstructMatch(ester_carbon_pattern) or 
            mol.HasSubstructMatch(amide_carbon_pattern)):
        return False, "Carbon connected to sulfur is not part of a long chain or ester/amide bond"

    return True, "Contains a sulfonic acid group connected to a lipid-like structure via a carbon-sulfur bond"