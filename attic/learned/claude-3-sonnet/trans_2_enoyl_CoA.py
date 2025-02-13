"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
            # Core CoA structure components
            "adenine": Chem.MolFromSmarts("c1nc(N)c2ncnc2n1"),
            "phosphate": Chem.MolFromSmarts("OP(=O)(O)O"),
            "thioester": Chem.MolFromSmarts("C(=O)SC"),
            "pantetheine": Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)"),
            # Pattern for alpha-beta unsaturated thioester, without stereochemistry
            "enoyl": Chem.MolFromSmarts("C=CC(=O)SC")
        }

        # Check for CoA moiety components
        if not mol.HasSubstructMatch(patterns["adenine"]):
            return False, "No CoA moiety found (missing adenine)"

        phosphate_matches = len(mol.GetSubstructMatches(patterns["phosphate"]))
        if phosphate_matches < 2:
            return False, "No CoA moiety found (insufficient phosphate groups)"

        if not mol.HasSubstructMatch(patterns["thioester"]):
            return False, "No thioester linkage found"

        if not mol.HasSubstructMatch(patterns["pantetheine"]):
            return False, "Missing pantetheine arm of CoA"

        # Find all double bonds connected to thioester
        enoyl_matches = mol.GetSubstructMatches(patterns["enoyl"])
        if not enoyl_matches:
            return False, "No α,β-unsaturated thioester found"

        # Check stereochemistry of the double bonds
        found_trans = False
        for match in enoyl_matches:
            # Get the double bond atoms
            double_bond_atoms = [match[0], match[1]]  # C=C atoms from the match
            bond = mol.GetBondBetweenAtoms(double_bond_atoms[0], double_bond_atoms[1])
            
            # Check if it's a double bond and has stereochemistry
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                # E/trans configuration is represented by STEREOCIS in RDKit
                if bond.GetStereo() == Chem.BondStereo.STEREOE:
                    found_trans = True
                    break

        if not found_trans:
            return False, "No trans configuration found at position 2"

        # Additional check for position (should be alpha to the thioester)
        # This is already enforced by our enoyl pattern

        return True, "Contains CoA moiety with trans-2-enoyl group"

    except Exception as e:
        return False, f"Error in structure analysis: {str(e)}"