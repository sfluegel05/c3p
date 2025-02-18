"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: CHEBI:64435 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol has a glycerol backbone with exactly two substituent groups
    (acyl, alkyl, or alk-1-enyl) attached via oxygen atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adjusted glycerol pattern to account for possible substituents
    # Three contiguous carbons with at least two CH2 groups (allowing any bonding)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # Define substituent patterns (acyl, alkyl, alk-1-enyl)
    ester = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")  # Acyl (ester)
    ether = Chem.MolFromSmarts("[OX2][CX4H2]")        # Alkyl (ether, at least two carbons)
    vinyl_ether = Chem.MolFromSmarts("[OX2]/C=C/")    # Alk-1-enyl (vinyl ether with trans configuration)
    vinyl_ether_cis = Chem.MolFromSmarts("[OX2]\\C=C\\")  # Cis vinyl ether

    # Check all possible glycerol backbones
    for backbone in matches:
        substituents = 0
        # Check each carbon in the backbone for substituents
        for carbon_idx in backbone:
            atom = mol.GetAtomWithIdx(carbon_idx)
            # Look for oxygen neighbors
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    # Check if oxygen is part of a substituent
                    o_bond = atom.GetBondBetweenAtoms(carbon_idx, neighbor.GetIdx())
                    if o_bond.GetBondType() != Chem.BondType.SINGLE:
                        continue  # Only single bonds for substituents
                    # Check for ester, ether, or vinyl ether patterns connected to this oxygen
                    context = Chem.FindAtomEnvironmentOfRadiusN(mol, radius=2, rootedAtAtomIdx=neighbor.GetIdx())
                    submol = Chem.PathToSubmol(mol, context)
                    if submol.HasSubstructMatch(ester):
                        substituents +=1
                        break  # One substituent per carbon
                    elif submol.HasSubstructMatch(ether) or submol.HasSubstructMatch(vinyl_ether) or submol.HasSubstructMatch(vinyl_ether_cis):
                        substituents +=1
                        break
        # Check for exactly two substituents
        if substituents == 2:
            return True, "Contains glycerol backbone with exactly two substituent groups (acyl/alkyl/alk-1-enyl)"
    
    return False, f"Found {substituents} substituent groups, need exactly 2"