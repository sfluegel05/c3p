"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def create_pattern(smarts):
    """Helper function to safely create and verify SMARTS patterns"""
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError(f"Invalid SMARTS pattern: {smarts}")
    return pattern

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    These are CoA thioesters with a trans double bond at position 2 of the acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Essential structural patterns
        patterns = {
            "adenine": "c1nc(N)c2ncnc2n1",
            "phosphate": "OP(=O)(O)O",
            "thioester": "C(=O)SC",
            "pantetheine": "SCCNC(=O)CCNC(=O)",
            # Two possible SMARTS representations for trans double bond
            "trans_enoyl_1": "\\C=C\\C(=O)S",
            "trans_enoyl_2": "/C=C/C(=O)S",
            "cis_enoyl_1": "/C=C\\C(=O)S",
            "cis_enoyl_2": "\\C=C/C(=O)S"
        }
        
        # Create all patterns
        mol_patterns = {name: create_pattern(smarts) for name, smarts in patterns.items()}

        # Check for CoA moiety components
        if not mol.HasSubstructMatch(mol_patterns["adenine"]):
            return False, "No CoA moiety found (missing adenine)"

        phosphate_matches = len(mol.GetSubstructMatches(mol_patterns["phosphate"]))
        if phosphate_matches < 2:
            return False, "No CoA moiety found (insufficient phosphate groups)"

        if not mol.HasSubstructMatch(mol_patterns["thioester"]):
            return False, "No thioester linkage found"

        if not mol.HasSubstructMatch(mol_patterns["pantetheine"]):
            return False, "Missing pantetheine arm of CoA"

        # Check for trans double bond configuration
        has_trans = (mol.HasSubstructMatch(mol_patterns["trans_enoyl_1"]) or 
                    mol.HasSubstructMatch(mol_patterns["trans_enoyl_2"]))
        has_cis = (mol.HasSubstructMatch(mol_patterns["cis_enoyl_1"]) or 
                  mol.HasSubstructMatch(mol_patterns["cis_enoyl_2"]))

        if not (has_trans or has_cis):
            return False, "No α,β-unsaturated thioester found"
        
        if has_cis and not has_trans:
            return False, "Found cis-enoyl-CoA instead of trans"

        # Additional checks for position 2
        # Count carbons between double bond and thioester to ensure position 2
        acyl_chain_pattern = create_pattern("C=CC(=O)S")
        if not mol.HasSubstructMatch(acyl_chain_pattern):
            return False, "Double bond not in position 2"

        return True, "Contains CoA moiety with trans-2-enoyl group"

    except Exception as e:
        return False, f"Error in structure analysis: {str(e)}"