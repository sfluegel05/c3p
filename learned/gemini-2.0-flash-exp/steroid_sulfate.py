"""
Classifies: CHEBI:16158 steroid sulfate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a steroid molecule with one or more sulfate groups attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general SMARTS pattern for the steroid core (tetracyclic ring system)
    # This pattern allows for a fused 6-6-6-5 ring system with various atom types in the rings. It is a generic representation.
    steroid_core_pattern = Chem.MolFromSmarts("[!H0;R3;R4;R5;R6]1[!H0;R3;R4;R5;R6][!H0;R3;R4;R5;R6]2[!H0;R3;R4;R5;R6]3[!H0;R3;R4;R5;R6]([!H0;R3;R4;R5;R6]1[!H0;R3;R4;R5;R6])[!H0;R3;R4;R5;R6][!H0;R3;R4;R5;R6]4[!H0;R3;R4;R5;R6]([!H0;R3;R4;R5;R6]2)([!H0;R3;R4;R5;R6][!H0;R3;R4;R5;R6]34)")
    if not mol.HasSubstructMatch(steroid_core_pattern):
         return False, "No steroid core found"
    
    # Define SMARTS pattern for the sulfate or sulfonate group, allowing for negative charges. 
    # The sulfate group is -O-S(=O)(=O)-O or -O-S(=O)(=O)-[O-]
    # Allow the sulfate group to be connected to any atom C, N, or O
    sulfate_pattern = Chem.MolFromSmarts("[#6,#7,#8]-[OX2]-[SX4](=[OX1])(=[OX1])([OX2,O-])")
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)

    if len(sulfate_matches) == 0:
        return False, "No sulfate group found"

    # Check if molecule contains at least 17 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 17:
        return False, "Steroid core must have at least 17 carbons."
    
     # Check molecular weight - steroid sulfates typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for steroid sulfate"

    return True, "Contains steroid core with one or more sulfate groups attached"