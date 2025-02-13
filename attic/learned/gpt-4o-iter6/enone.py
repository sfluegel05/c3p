"""
Classifies: CHEBI:51689 enone
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is an alpha, beta-unsaturated ketone with the C=O function
    conjugated to a C=C double bond at the alpha, beta position, and R(4) not being hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern to match enone structures: C=O conjugated to C=C
    # The pattern ensures the C=C bond in conjugation with the C=O group
    enone_pattern = Chem.MolFromSmarts("C=CC=O")

    # Check if the molecule has the enone substructure
    if mol.HasSubstructMatch(enone_pattern):
        # Further verify that the carbonyl carbon is not attached to hydrogen
        # R1RC=CR-C(=O)R4, where R4 != H
        carbonyl_oxygen = "[CX3]=[OX1]"    # Standard carbonyl group
        ketone_check_pattern = Chem.MolFromSmarts(carbonyl_oxygen)
        match = mol.GetSubstructMatch(ketone_check_pattern)

        if match:
            carbon = match[0]   # Get the carbonyl carbon
            neighbor_atoms = mol.GetAtomWithIdx(carbon).GetNeighbors()

            # Ensure none of the neighbors is hydrogen
            for atom in neighbor_atoms:
                if atom.GetAtomicNum() == 1:
                    return False, "Carbonyl carbon in SMILES structure is bonded to hydrogen"

            return True, "Contains alpha, beta-unsaturated ketone (enone) structure"
            
    return False, "No enone structure found"