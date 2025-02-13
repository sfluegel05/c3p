"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: CHEBI:26194 sphingoid
Definition: Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sphinganine core (long aliphatic chain with amine group)
    core_pattern = Chem.MolFromSmarts("[CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CX4]([CH2X4][CH2X4][CH2X4][CH2X4][CH3X4])[NX3H2,NX3,NX4H2+]")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No sphinganine core found"
    
    # Look for hydroxy groups
    has_hydroxy = any(atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1 for atom in mol.GetAtoms())
    
    # Look for unsaturated bonds
    has_unsaturation = any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in mol.GetBonds())
    
    # Sphingoids must have either hydroxy groups or unsaturated bonds
    if not (has_hydroxy or has_unsaturation):
        return False, "No hydroxy groups or unsaturated bonds found"
    
    # Check for chiral centers (optional)
    chiral_centers = any(atom.HasProp("_ChiralityPossible") and atom.HasProp("_CIPCode") and atom.GetProp("_CIPCode") in ["R", "S"] for atom in mol.GetAtoms())
    
    return True, "Contains sphinganine core with hydroxy/unsaturated groups" + (", chiral centers present" if chiral_centers else "")