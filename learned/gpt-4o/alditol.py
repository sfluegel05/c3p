"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with the general formula HOCH2-[CH(OH)]n-CH2OH, 
    derived from an aldose by reduction of the carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alditol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an alditol structure
    # Correct pattern for HOCH2-[CH(OH)]n-CH2OH might look like "OCC(O)C(O)C(O)*CO"
    # However, we need to ensure we are matching linear structures, hence a better
    # pattern would target all hydroxylated carbon chains of sufficient length
    alditol_pattern = Chem.MolFromSmarts("[CX4](O)CO")  # Generic repeating unit for acyclic polyols
    if alditol_pattern is None:
        return (None, "Invalid SMARTS pattern")

    # Check if the molecule contains cycles
    if rdmolops.GetSSSR(mol) > 0:
        return False, "Molecule contains a ring structure, not an acyclic alditol"

    # Calculate the total number of carbon atoms
    carbon_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    # Check if molecule has enough contiguous polyol chain as per alditol requirements
    # Minimum of 3 pairs of [C-O] for HO(CH(OH))nCH2OH with n >= 2 and ends with -CH2OH
    match = mol.GetSubstructMatches(alditol_pattern)
    if len(match) < 2 or carbon_atoms < 4:
        return False, f"Structure with {len(match)} C-O patterns, insufficient for alditol"
    
    return True, "SMILES string matches the structural pattern of an acyclic alditol"