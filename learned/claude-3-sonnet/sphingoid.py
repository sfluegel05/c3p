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
    
    # Look for sphinganine backbone pattern (long carbon chain with amine and 2+ hydroxy groups)
    backbone_pattern = Chem.MolFromSmarts("[CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CH2X4][CX4]([CH2X4][CH2X4][CH2X4][CH2X4][CH3X4])[NX3H2,NX3,NX4H2+][CX4]([OX2H,OX1H,OX1-])[CX4]([OX2H,OX1H,OX1-])")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No sphinganine backbone found"
    
    # Look for long aliphatic chains
    chain_pattern = Chem.MolFromSmarts("[CH2X4,CH3X4]~[CH2X4,CH3X4]~[CH2X4,CH3X4]~[CH2X4,CH3X4]~[CH2X4,CH3X4]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 2:
        return False, "Missing long aliphatic chains"
    
    # Check for hydroxy groups or unsaturated bonds
    has_hydroxy = any(atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1 for atom in mol.GetAtoms())
    has_unsaturation = any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in mol.GetBonds())
    if not (has_hydroxy or has_unsaturation):
        return False, "No hydroxy groups or unsaturated bonds found"
    
    # Check for chiral centers (optional)
    chiral_centers = any(atom.HasProp("_ChiralityPossible") and atom.HasProp("_CIPCode") and atom.GetProp("_CIPCode") in ["R", "S"] for atom in mol.GetAtoms())
    
    return True, "Contains sphinganine backbone with long aliphatic chains and hydroxy/unsaturated groups" + (", chiral centers present" if chiral_centers else "")