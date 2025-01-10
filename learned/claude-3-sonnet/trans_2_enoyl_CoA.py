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
            # Core CoA structure with complete specificity
            "coA_core": Chem.MolFromSmarts("[$(N1C=NC2=C(N)N=CN=C12)]"),
            
            # Specific pattern for trans-2-enoyl-CoA:
            # - Requires thioester (SC(=O))
            # - Connected to exactly one trans double bond at position 2
            # - R group can be any carbon chain
            "trans_2_enoyl": Chem.MolFromSmarts("[#6]/[CH]=[CH]/C(=O)SC"),
            
            # Pattern to exclude compounds with additional conjugated double bonds
            "conjugated_diene": Chem.MolFromSmarts("[#6]/[CH]=[CH]/[CH]=[CH]/[#6]"),
        }

        # Check for CoA core
        if not mol.HasSubstructMatch(patterns["coA_core"]):
            return False, "Missing CoA core structure"

        # Find trans-2-enoyl pattern matches
        trans_2_enoyl_matches = mol.GetSubstructMatches(patterns["trans_2_enoyl"])
        
        if not trans_2_enoyl_matches:
            return False, "No trans-2-enoyl group found"

        # Check for conjugated dienes (should not be present)
        if mol.HasSubstructMatch(patterns["conjugated_diene"]):
            # Get the positions of conjugated double bonds
            diene_matches = mol.GetSubstructMatches(patterns["conjugated_diene"])
            # Check if any of these overlap with our trans-2-enoyl position
            for diene_match in diene_matches:
                for trans_match in trans_2_enoyl_matches:
                    if set(diene_match[:4]).intersection(set(trans_match[:2])):
                        return False, "Contains conjugated diene system"

        # Verify double bond stereochemistry
        for match in trans_2_enoyl_matches:
            # Get the double bond atoms
            double_bond_atoms = [match[1], match[2]]  # C=C atoms from the match
            bond = mol.GetBondBetweenAtoms(double_bond_atoms[0], double_bond_atoms[1])
            
            # Must be a double bond with explicit trans (E) stereochemistry
            if (bond.GetBondType() != Chem.BondType.DOUBLE or 
                bond.GetStereo() != Chem.BondStereo.STEREOE):
                return False, "Double bond must have explicit trans (E) stereochemistry"

            # Check that this is the only double bond connected to the thioester
            thioester_carbon = match[3]  # C(=O) carbon
            for neighbor in mol.GetAtomWithIdx(thioester_carbon).GetNeighbors():
                if neighbor.GetIdx() not in [match[2], match[4]]:  # not the C= or S atoms
                    return False, "Additional substitution at thioester position"

        return True, "Contains CoA moiety with trans-2-enoyl group"

    except Exception as e:
        return False, f"Error in structure analysis: {str(e)}"