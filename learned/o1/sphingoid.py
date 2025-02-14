"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: sphingoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid is defined as sphinganine, its homologs and stereoisomers,
    and the hydroxy and unsaturated derivatives of these compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sphingoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for key features

    # Long aliphatic chain (at least 8 contiguous carbons)
    chain_pattern = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long aliphatic chain found (minimum 8 carbons)"

    # Amino alcohol backbone
    # Nitrogen atom connected to a carbon with at least one hydroxyl group
    amino_alcohol_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(N[C,S]=O)]-[CX4]-[OX2H]")
    if not mol.HasSubstructMatch(amino_alcohol_pattern):
        return False, "No amino alcohol backbone found"

    # Additional hydroxyl group (may be at C3 position)
    # Oxygen atom connected to a carbon adjacent to the amino group
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]-[CX4]-[CX4]-[NX3;H2,H1;!$(N[C,S]=O)]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No additional hydroxyl group found"

    # Check for possible unsaturation (double bonds) or hydroxyl derivatives
    # This is allowed, so we don't need to explicitly check for them

    # Count total number of carbons (should be at least 12)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, f"Too few carbons ({c_count} found, minimum 12 required)"

    # Ensure the molecule is a long-chain amino diol
    # Find the path from amino group to terminal hydroxyl group
    # This ensures the backbone is continuous

    # Identify amino group nitrogen
    amino_n = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetTotalDegree() == 3]
    if not amino_n:
        return False, "No appropriate amino group found"
    amino_n = amino_n[0]

    # Identify primary alcohol oxygens
    primary_oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1 and atom.GetTotalNumHs() == 1]
    if len(primary_oxygens) < 1:
        return False, "No primary hydroxyl group found"

    # Check connectivity between amino group and hydroxyl groups
    # Use a BFS to find paths
    found_backbone = False
    for oxygen in primary_oxygens:
        # Get the carbon attached to the hydroxyl group
        oxy_carbon = oxygen.GetNeighbors()[0]
        # Find path to amino nitrogen
        path = Chem.rdmolops.GetShortestPath(mol, amino_n.GetIdx(), oxy_carbon.GetIdx())
        # Backbone should be at least 2 atoms (C-N)
        if len(path) >= 2:
            found_backbone = True
            break

    if not found_backbone:
        return False, "No proper amino alcohol backbone found"

    return True, "Molecule matches sphingoid structural features"