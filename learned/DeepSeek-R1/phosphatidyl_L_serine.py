"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:16015 phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    A phosphatidyl-L-serine consists of a glycerol backbone with two fatty acid esters,
    a phosphate group, and an L-serine moiety attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for phosphate group (O=P(O)(O)O...)
    phosphate_pattern = Chem.MolFromSmarts("[O][P](=O)([O])[O]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group"

    # Check for serine moiety: NH2-CH(COO)-OPO3 connected in proper structure
    # L-serine has S configuration at the alpha carbon (C4 in IUPAC numbering)
    serine_pattern = Chem.MolFromSmarts("[NH2,NH3+][C@](C(=O)[O-,OH])([OH,O-])OP(=O)([O-])[O-]")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "Missing L-serine moiety with correct stereochemistry"

    # Check for two ester groups (fatty acids)
    ester_pattern = Chem.MolFromSmarts("[OX2]C(=O)")
    esters = mol.GetSubstructMatches(ester_pattern)
    if len(esters) < 2:
        return False, f"Expected at least 2 ester groups, found {len(esters)}"

    # Check for glycerol backbone (O-C-O-C-O-C structure)
    glycerol_backbone = Chem.MolFromSmarts("[CH2]O[CH2]")
    if not mol.HasSubstructMatch(glycerol_backbone):
        return False, "Glycerol backbone not found"

    # Optional: Check approximate chain lengths (at least 8 carbons in each fatty acid)
    # This is a heuristic and may need adjustment
    long_chain = Chem.MolFromSmarts("C-C-C-C-C-C-C-C")
    if len(mol.GetSubstructMatches(long_chain)) < 2:
        return False, "Insufficient long-chain fatty acids"

    return True, "Phosphatidyl group esterified to L-serine with two fatty acid chains"