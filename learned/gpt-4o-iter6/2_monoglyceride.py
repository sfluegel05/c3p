"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride is defined by an ester linkage at the second carbon of a glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 2-monoglyceride pattern:
    # [CX4] = central carbon with 4 attachments
    # Ester linkage (O-C(=O)) at the second carbon of glycerol [-CH(OC(=O)R)-]
    # Hydroxyls (-OH) at the first and third positions [-CHOH-]
    glycerol_2_mono_pattern = Chem.MolFromSmarts("OCC(O)C(OC(=O))")
    
    if not mol.HasSubstructMatch(glycerol_2_mono_pattern):
        return False, "No specific 2-monoglyceride pattern found with esterification at C2"

    return True, "Contains 2-monoglyceride structure with acyl group esterified at the second position"