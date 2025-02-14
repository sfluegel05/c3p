"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: CHEBI:36838 Fatty acid anion
The conjugate base of a fatty acid, arising from deprotonation of the carboxylic acid group of the corresponding fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylate group (-COO-)
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Look for long carbon chain (>=12 carbons)
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No long carbon chain found"

    # Disqualifying functional groups
    disqualifying_patterns = [
        Chem.MolFromSmarts("[#6]=[#6]"),  # Carbon-carbon double bonds
        Chem.MolFromSmarts("[#6]#[#6]"),  # Carbon-carbon triple bonds
        Chem.MolFromSmarts("[!#6;!#1]"),  # Non-carbon, non-hydrogen atoms
        Chem.MolFromSmarts("[#8]~[#8]"),  # Peroxides
        Chem.MolFromSmarts("[#6]-[#6](-[#8])-[#6]"),  # Ethers
        Chem.MolFromSmarts("[#6]-[#6](-[#6])-[#6]"),  # Tertiary carbons
        Chem.MolFromSmarts("[#6]=,:[#6]"),  # Aromatic rings
    ]

    for pattern in disqualifying_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, f"Contains disqualifying functional group: {pattern.GetSmarts()}"

    # Additional structural constraints
    # 1. Carboxylate group must be at the end of the carbon chain
    chain_atoms = [atom for match in chain_matches for atom in mol.GetAtomWithIdx(match)]
    if not any(atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1 for atom in chain_atoms):
        return False, "Carboxylate group not at the end of the carbon chain"

    # 2. No rings or cyclic structures
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains ring or cyclic structure"

    # 3. No branched chains
    for atom in chain_atoms:
        if atom.GetDegree() > 2:
            return False, "Contains branched carbon chain"

    # Consider stereochemistry
    # Stereochemistry is important for some fatty acid anions
    # For simplicity, we will not consider stereochemistry in this implementation

    # Handle exceptional cases
    # There may be some exceptional cases or edge cases that the program is not handling correctly
    # We will address these as they are encountered

    # Check molecular weight - typical range for fatty acid anions
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, "Molecular weight outside typical range for fatty acid anions"

    return True, "Molecule meets the structural requirements for a fatty acid anion"