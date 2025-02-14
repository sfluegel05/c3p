"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a tetrasaccharide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for pyranose and furanose rings 
    # including stereochemistry-neutral patterns
    pyranose_pattern = Chem.MolFromSmarts("[C&O][C&OH]1O[C&OH][C&OH][C&OH][C&OH]1")
    furanose_pattern = Chem.MolFromSmarts("[C&O][C&OH]1OC[C&OH][C&OH]1")
    
    # Count matches for pyranose and furanose rings
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)

    # Each ring implies a monosaccharide unit
    num_monosaccharide_units = len(pyranose_matches) + len(furanose_matches)
    
    if num_monosaccharide_units == 4:
        return True, "Contains four linked monosaccharide units"
    else:
        return False, f"Found {num_monosaccharide_units} monosaccharide units, need exactly 4"

# Example usage - Replace with actual SMILES for testing
# print(is_tetrasaccharide("some_smiles_string"))