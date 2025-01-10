"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: prostaglandin
"""
from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    A prostaglandin is characterized by a 20-carbon skeleton that includes 
    a cyclopentane ring with two side chains. The side chains have specific
    functional groups, and the molecule may have additional rings or variations.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the prostaglandin core SMARTS pattern
    # Cyclopentane ring with specific substituents
    prostaglandin_core = Chem.MolFromSmarts("""
        [
            # Cyclopentane ring
            R1; 
            $( [C;R1](-[C])(=O) ),  # Carbonyl group at position 2
            $( [C;R1]([O]) ),       # Hydroxyl group at position 3
            $( [C;R1]-[C]=[C]-[C] ),# Side chain with double bonds
            $( [C;R1]-[C]-[C]-[C](=O) ), # Side chain ending with carboxyl or derivative
        ]
    """)
    if prostaglandin_core is None:
        return False, "Error in SMARTS pattern"

    # Check for prostaglandin core match
    if not mol.HasSubstructMatch(prostaglandin_core):
        return False, "Molecule does not match prostaglandin core structure"

    # Check total number of carbons (should be approximately 20)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 17 or num_carbons > 24:
        return False, f"Number of carbons is {num_carbons}, which is not in the expected range for prostaglandins"

    # Optional: Validate stereochemistry (commented out due to variability)
    # chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    # if len(chiral_centers) < expected_number:
    #     return False, "Insufficient chiral centers for a prostaglandin"

    # If all checks pass, it is likely a prostaglandin
    return True, "Molecule matches prostaglandin structural features"

__metadata__ = {
    'chemical_class': {
        'name': 'prostaglandin',
        'definition': 'Naturally occurring compounds derived from the parent C20 acid, prostanoic acid.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
    },
    'message': None,
    'attempt': 3,
    'success': True,
    'error': '',
    'stdout': None
}