"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    A mineral is generally a chemical substance that is normally crystalline and formed through geological processes.
    This includes inorganic compounds, salts, and certain naturally occurring amorphous substances.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Common mineral-forming metals (including alkali, alkaline earth, transition metals, etc.)
    metal_atomic_nums = set(range(3, 32)) | {37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
                                            55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68,
                                            69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82,
                                            83, 88, 89, 90, 91, 92}

    # Check for the presence of metal ions (common in minerals)
    metal_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in metal_atomic_nums]
    if not metal_atoms:
        return False, "No metal ions found, which are common in minerals"

    # Check for common mineral patterns with more flexible matching
    patterns = {
        "sulfate": "[O-]S(=O)(=O)[O-]",
        "carbonate": "[O-]C(=O)[O-]",
        "phosphate": "[O-]P(=O)([O-])[O-]",
        "silicate": "[O-][Si]([O-])([O-])[O-]",
        "oxide": "[O-2]",
        "sulfide": "[S-2]",
        "halide": "[F,Cl,Br,I]-",
        "hydroxide": "[OH-]",
        "nitrate": "[O-][N+]([O-])=O",
        "water": "O"
    }

    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return True, f"Contains {name} group, common in minerals"

    # Check for simple inorganic salts
    if len(mol.GetAtoms()) <= 5 and all(atom.GetAtomicNum() in metal_atomic_nums or 
                                      atom.GetAtomicNum() in {6, 7, 8, 9, 15, 16, 17, 35, 53} for atom in mol.GetAtoms()):
        return True, "Simple inorganic salt, common in minerals"

    # If none of the above patterns match, return False
    return False, "Does not match common mineral patterns"